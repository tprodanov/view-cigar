#!/usr/bin/env python3

from collections import namedtuple
import sys
import builtins
import re
import logging
import pysam
from colored import Fore, Back, Style


def panic(message):
    logging.critical(message)
    exit(1)


def open(filename, mode):
    if filename == '-':
        return sys.stdin if 'r' in mode else sys.stdout
    elif filename.endswith('.gz'):
        return gzip.open(filename, mode)
    elif filename.endswith('.xz'):
        return lzma.open(filename, mode)
    return builtins.open(filename, mode)


Alignment = namedtuple('Alignment',
    ['ref_name', 'ref_start', 'ref_end',
    'query_name', 'query_start', 'query_end', 'query_strand',
    'cigar', 'query_seq'])


def load_paf(f):
    for line in f:
        if line.startswith('#') or line.startswith('@'):
            continue
        parts = line.strip().split('\t')
        query_name = parts[0]
        query_start = int(parts[2])
        query_end = int(parts[3])
        query_strand = parts[4] == '+'
        ref_name = parts[5]
        ref_start = int(parts[7])
        ref_end = int(parts[8])
        cigar = None
        for entry in parts[12:]:
            if entry.startswith('cg:Z:'):
                cigar = entry[5:]
                break
        if cigar is None:
            logging.warn(f'No CIGAR found for {query_name} and {ref_name}, skipping')
        else:
            yield Alignment(ref_name, ref_start, ref_end, query_name, query_start, query_end, query_strand, cigar, None)


def load_alignments(args):
    assert args.cigar is not None or args.paf is not None or args.sam is not None
    if args.cigar is not None:
        if 'H' in args.cigar:
            panic('Hard clipping in CIGAR is only allowed for SAM files')
        yield Alignment('Reference', 0, None, 'Query', 0, None, True, args.cigar, None)

    elif args.paf is not None:
        with open(args.paf, 'rt') as f:
            yield from load_paf(f)

    elif args.sam is not None:
        raise NotImplementedError()


class ColoredString:
    def __init__(self):
        self.text = ''
        self.active_code = ''

    def extend(self, text, code=''):
        if self.active_code != code:
            if self.active_code != '':
                self.text += Style.reset
            self.text += code
            self.active_code = code
        self.text += str(text)

    def __iadd__(self, text):
        # Add text without changing the code.
        self.text += str(text)
        return self

    def finish(self):
        if self.active_code != '':
            self.text += Style.reset
            self.active_code = ''

    def __str__(self):
        self.finish()
        return self.text


class PlainString:
    def __init__(self):
        self.text = ''

    def extend(self, text, _code=''):
        self.text += str(text)

    def __iadd__(self, text):
        self.text += str(text)
        return self

    def finish(self):
        pass

    def __str__(self):
        return self.text


class Painter:
    def __init__(self, color_mode, stream):
        if color_mode == 'always':
            self.colored = True
        elif color_mode == 'never':
            self.colored = False
        else:
            self.colored = stream.isatty()

    def __call__(self, text, *codes):
        if self.colored and codes:
            return f'{"".join(codes)}{text}{Style.reset}'
        return str(text)

    def new_string(self):
        return ColoredString() if self.colored else PlainString()


class Sequences:
    def __init__(self):
        self.indexed_fasta = None
        self.sequences = None
        self.sequence = None

    def fetch(self, name, start, end):
        if self.indexed_fasta is not None:
            try:
                return self.indexed_fasta.fetch(name, start, end).upper()
            except (IndexError, ValueError) as e:
                logging.error(f'Cannot fetch {name}:{start+1}-{end} (sequence not found or incorrect coordinates)')
                return None

        seq = self.sequence if self.sequence else self.sequences.get(name)
        if seq is None:
            logging.error(f'Cannot find sequence {name}')
            return None
        if end is None:
            return seq
        if len(seq) < end:
            logging.error(f'Sequence {name} is too short ({len(seq)}, '
                f'asked for region {start+1}-{end})')
            return None
        return seq[start:end]


def load_sequences(filename, sequence, tag):
    if filename is not None and sequence is not None:
        panic(f'Both {tag} file and sequence were provided')
    if filename is None and sequence is None:
        panic(f'Neither {tag} file nor sequence were provided')

    res = Sequences()
    if sequence is not None:
        res.sequence = sequence.upper()
        return res

    try:
        fasta = pysam.FastaFile(filename)
        res.indexed_fasta = fasta
    except ValueError:
        # No index found.
        res.sequences = {}
        with pysam.FastxFile(filename) in fastx:
            for record in fastx:
                res.sequences[record.name] = record.sequence.upper()
    return res


OPERATION_REGEX = re.compile(r'([0-9]+)([MIDNSHP=X])')

def process_cigar(cigar):
    """
    Processes CIGAR and returns a list of pairs (query_pos, ref_pos),
    as well as the alignment query and reference lengths.
    """
    ref_pos = 0
    query_pos = 0
    pairs = []
    for length, operation in OPERATION_REGEX.findall(cigar):
        length = int(length)
        match operation:
            case 'M' | '=' | 'X':
                for i in range(length):
                    pairs.append((query_pos + i, ref_pos + i))
                ref_pos += length
                query_pos += length
            case 'I' | 'S' | 'H':
                # Put hard padding here because we will replace cropped read sequence with full sequence
                for i in range(length):
                    pairs.append((query_pos + i, None))
                query_pos += length
            case 'D' | 'N':
                for i in range(length):
                    pairs.append((None, ref_pos + i))
                ref_pos += length
            case 'P':
                for i in range(length):
                    pairs.append((None, None))
    return pairs, query_pos, ref_pos


LETTER_COLORS = { 'A': Fore.green, 'C': Fore.blue, 'G': Fore.yellow, 'T': Fore.red, '-': Fore.light_gray }


def visualize_alignment(aln, ref_seqs, query_seqs, out, painter, width):
    ref_seq = ref_seqs.fetch(aln.ref_name, aln.ref_start, aln.ref_end)
    query_seq = aln.query_seq if aln.query_seq is not None \
        else query_seqs.fetch(aln.query_name, aln.query_start, aln.query_end)
    if ref_seq is None or query_seq is None:
        out.write(painter(f'\n\n!!! Unavailable sequences "{aln.query_name}" and "{aln.ref_name}" !!!\n\n', Fore.red))
        return
    if not aln.query_strand:
        query_seq = pysam.reverse_complement(query_seq)

    pairs, query_len, ref_len = process_cigar(aln.cigar)
    if len(ref_seq) != ref_len or len(query_seq) != query_len:
        logging.error(f'Alignment {aln.query_name} to {aln.ref_name}: incorrect alignment length: '
            f'reference {ref_len} (expected {len(ref_seq)}), query {query_len} (expected {len(query_seq)})')
        out.write(painter(f'\n\n!!! Incorrect alignment {aln.query_name} to {aln.ref_name} !!!\n\n', Fore.red))
        return

    out.write('=====\n')
    out.write(painter(f'{aln.query_name}:{aln.query_start+1:,}-{aln.query_start + query_len:,}', Style.bold))
    out.write(' (top) to ')
    out.write(painter(f'{aln.ref_name}:{aln.ref_start+1:,}-{aln.ref_start + ref_len:,}', Style.bold))
    out.write(' (bottom)\n')

    max_end_length = len('{:,}'.format(max(aln.ref_start + ref_len, aln.query_start + query_len)))
    left_margin = max_end_length - 1
    margin = ' ' * left_margin
    if left_margin > width + 10:
        width = left_margin + 10
        logging.warning(f'Width too small, expanding to {width}')
    active_width = width - left_margin

    lines = []
    offset = 0
    query_start = aln.query_start if aln.query_strand else -(aln.query_end - 1)
    ref_start = aln.ref_start
    last_qpos = abs(query_start) + 1
    last_rpos = ref_start

    for i in range(0, len(pairs), active_width):
        margin_q_offset, margin_r_offset = pairs[i]
        margin_qpos = abs(query_start + margin_q_offset) + 1 if margin_q_offset is not None else last_qpos
        margin_rpos = ref_start + margin_r_offset + 1 if margin_r_offset is not None else last_rpos
        qpos_line = f'{margin_qpos:>{max_end_length},}\n'
        rpos_line = f'{margin_rpos:>{max_end_length},}\n'
        qdigit_line = painter.new_string()
        rdigit_line = painter.new_string()
        aln_line = ''
        qletter_line = painter.new_string()
        rletter_line = painter.new_string()

        for query_offset, ref_offset in pairs[i:i + active_width]:
            if query_offset is None:
                qdigit_line += ' '
                qletter = '-'
            else:
                last_qpos = abs(query_start + query_offset) + 1
                digit = last_qpos % 10
                qdigit_line.extend(digit, '' if digit == 0 else Fore.light_gray)
                qletter = query_seq[query_offset]

            if ref_offset is None:
                rdigit_line += ' '
                rletter = '-'
            else:
                last_rpos = ref_start + ref_offset + 1
                digit = last_rpos % 10
                rdigit_line.extend(digit, '' if digit == 0 else Fore.light_gray)
                rletter = ref_seq[ref_offset]

            aln_line += ' ' if qletter != rletter or query_offset is None else '|'
            qletter_line.extend(qletter, LETTER_COLORS.get(qletter, ''))
            rletter_line.extend(rletter, LETTER_COLORS.get(rletter, ''))

        out.write('\n')
        out.write(qpos_line)
        out.write(f'{margin}{qdigit_line}\n')
        out.write(f'{margin}{qletter_line}\n')
        out.write(f'{margin}{aln_line}\n')
        out.write(f'{margin}{rletter_line}\n')
        out.write(f'{margin}{rdigit_line}\n')
        out.write(rpos_line)


def main():
    import argparse
    import importlib.metadata
    parser = argparse.ArgumentParser(
        description='Visualize alignment.',
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False,
        usage='%(prog)s (-s SAM | -p PAF | -c CIGAR) (-R FASTA -Q FASTA | -r SEQ -q SEQ) [options]')
    aln_args = parser.add_argument_group('Alignment')
    aln_group = aln_args.add_mutually_exclusive_group(required=True)
    aln_group.add_argument('-s', '--sam', metavar='FILE',
        help='Input SAM/BAM file. Mutually exclusive with -p and -c.\n'
            'Does not require query sequences (-Q/-q).')
    aln_group.add_argument('-p', '--paf', metavar='FILE',
        help='Input PAF file. Mutually exclusive with -s and -c.')
    aln_group.add_argument('-c', '--cigar', metavar='STR',
        help='Cigar string. Mutually exclusive with -s and -p.\n'
            'Must match reference and query sequences from start to end\n'
            'on a positive strand.')

    seq_args = parser.add_argument_group('Sequences')
    seq_args.add_argument('-R', '--ref-fa', metavar='FILE',
        help='Optionally indexed FASTA/Q file. Must contain reference sequences.')
    seq_args.add_argument('-Q', '--query-fa', metavar='FILE',
        help='Optionally indexed FASTA/Q file. Must contain query sequences.\n'
            'Equals to -R if -Q/-q are not provided.')
    seq_args.add_argument('-r', '--ref-seq', metavar='STR',
        help='Reference sequence. Mutually exclusive with -f, requires -q.')
    seq_args.add_argument('-q', '--query-seq', metavar='STR',
        help='Reference sequence. Mutually exclusive with -f, requires -r.')

    opt_args = parser.add_argument_group('Optional arguments')
    opt_args.add_argument('-o', '--output', metavar='FILE', default='-',
        help='Output file [stdout].')
    opt_args.add_argument('-w', '--width', metavar='INT', default=120, type=int,
        help='Display width [%(default)s].')
    opt_args.add_argument('--color', metavar='WHEN', default='auto',
        help='When to provide colored output (always/auto/never) [%(default)s].')
    opt_args.add_argument('-h', '--help', action='help',
        help='Show this help message and exit.')
    opt_args.add_argument('-V', '--version', action='version',
        help="Show program's version number and exit.",
        version='visualize-cigar {version}'.format(version=importlib.metadata.version('visualize-cigar')))
    args = parser.parse_args()

    logger = logging.basicConfig(format='%(levelname)-7s| %(message)s', level=logging.DEBUG)

    ref_seqs = load_sequences(args.ref_fa, args.ref_seq, 'reference')
    if args.sam is not None:
        query_seqs = None
    elif args.query_fa is not None or args.query_seq is not None:
        query_seqs = load_sequences(args.query_fa, args.query_seq, 'query')
    elif args.ref_fa is not None:
        query_seqs = ref_seqs
    else:
        panic('At least one of -R/-Q/-q arguments must be provided.')

    out = open(args.output, 'wt')
    painter = Painter(args.color, out)
    for aln in load_alignments(args):
        visualize_alignment(aln, ref_seqs, query_seqs, out, painter, args.width)


if __name__ == '__main__':
    main()
