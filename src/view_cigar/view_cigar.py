#!/usr/bin/env python3

from collections import namedtuple, deque
import sys
import builtins
import re
import os, os.path
import contextlib
import shutil
import subprocess
import logging
import pysam
from colored import Fore, Back, Style


try:
    from pysam import reverse_complement
except ImportError:
    REV_COMP_NT = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A' }
    def reverse_complement(seq):
        return ''.join(REV_COMP_NT.get(nt, 'N') for nt in seq[::-1])


@contextlib.contextmanager
def paginator(enabled):
    cmd = ['less', '-SRFX']
    """Yield a file-like object that is piped through a pager if enabled."""
    if not enabled or not sys.stdout.isatty() or shutil.which(cmd[0]) is None:
        yield sys.stdout
        return

    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)
    try:
        proc.stdin.isatty = lambda: True
        yield proc.stdin
    except (BrokenPipeError, KeyboardInterrupt) as e:
        pass
    finally:
        try:
            proc.stdin.close()
        except BrokenPipeError:
            pass
        proc.wait()


def panic(message):
    logging.critical(message)
    exit(1)


def open(filename, mode, paginate=False):
    if filename == '-':
        return sys.stdin if 'r' in mode else paginator(paginate)
    elif filename.endswith('.gz'):
        return gzip.open(filename, mode)
    elif filename.endswith('.xz'):
        return lzma.open(filename, mode)
    return builtins.open(filename, mode)


class Interval:
    def __init__(self, name, start, end):
        self.name = name
        self.start = start
        self.end = end

    @classmethod
    def parse(cls, s):
        name, coords = s.split(':')
        start, end = map(int, coords.split('-'))
        return cls(name, start - 1, end)


class Alignment(namedtuple('Alignment', ['ref_interval', 'query_interval', 'query_strand', 'cigar', 'query_seq'])):
    @property
    def ref_name(self):
        return self.ref_interval.name

    @property
    def ref_start(self):
        return self.ref_interval.start

    @property
    def ref_end(self):
        return self.ref_interval.end

    @property
    def query_name(self):
        return self.query_interval.name

    @property
    def query_start(self):
        return self.query_interval.start

    @property
    def query_end(self):
        return self.query_interval.end


def load_sam(filename, ref_fasta):
    with pysam.AlignmentFile(filename, require_index=False,
            reference_filename=ref_fasta if filename.endswith('.cram') else None) as sam:
        last_name = None
        last_re = None
        last_seq = None
        last_strand = True
        for record in sam:
            name = record.query_name
            read_end = record.is_read1
            strand = record.is_forward
            cigar = record.cigarstring

            seq = record.query_sequence
            if seq is not None and 'H' not in cigar:
                last_name = name
                last_re = read_end
                last_strand = strand
                last_seq = seq
            else:
                if name != last_name or last_re != read_end:
                    logging.debug(f'Cannot identify sequence for read {name}, skipping')
                    continue
                seq = last_seq if strand == last_strand else reverse_complement(last_seq)
                cigar = cigar.replace('H', 'S')
            # Set strand to True because we already reverse-complemented sequence, if needed.
            yield Alignment(
                Interval(record.reference_name, record.reference_start, record.reference_end),
                Interval(name, 0, len(seq)), True, cigar, seq)


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
        elif 'H' in cigar:
            logging.error('Hard clipping in CIGAR is only allowed for SAM files, skipping')
        else:
            yield Alignment(
                Interval(ref_name, ref_start, ref_end),
                Interval(query_name, query_start, query_end),
                query_strand, cigar, None)


def load_alignments(args, force_ref_region, force_query_region, force_query_strand):
    assert args.cigar is not None or args.paf is not None or args.sam is not None
    if args.cigar is not None:
        if 'H' in args.cigar:
            logging.error('Hard clipping in CIGAR is only allowed for SAM files, skipping')
        else:
            if force_ref_region:
                ref_interval = force_ref_region
            else:
                ref_interval = Interval('Reference', 0, None)
            if force_query_region:
                query_interval = force_query_region
            else:
                query_interval = Interval('Query', 0, None)
            query_strand = True if force_query_strand is None else force_query_strand
            yield Alignment(ref_interval, query_interval, query_strand, args.cigar, None)

    elif args.paf is not None:
        with open(args.paf, 'rt') as f:
            yield from load_paf(f)

    elif args.sam is not None:
        yield from load_sam(args.sam, args.ref_fa)


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
        match color_mode.lower():
            case 'always' | 'yes' | 'y':
                self.colored = True
            case 'never' | 'n':
                self.colored = False
            case 'auto':
                self.colored = stream.isatty()
            case _:
                panic(f'Unsupported option --color {color_mode}')

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

    def fetch(self, name, start, end, padding):
        """
        Tries to get name:start-end, optionally extending ends by at most `padding`.
        Returns left_padding, right_padding, fetched_sequence.
        """
        if self.indexed_fasta is not None:
            try:
                seq_len = self.indexed_fasta.get_reference_length(name)
            except KeyError:
                logging.error(f'Cannot fetch {name}:{start+1}-{end} (sequence not found)')
                return None
            if seq_len < end:
                logging.error(f'Sequence {name} is too short ({seq_len}), '
                    f'asked for region {start+1}-{end})')
                return None
            left_padding = min(padding, start)
            right_padding = min(padding, seq_len - end)
            seq = self.indexed_fasta.fetch(name, start - left_padding, end + right_padding).upper()
            return left_padding, right_padding, seq

        seq = self.sequence if self.sequence else self.sequences.get(name)
        if seq is None:
            logging.error(f'Cannot find sequence {name}')
            return None
        if end is None:
            end = len(seq)
        if len(seq) < end:
            logging.error(f'Sequence {name} is too short ({len(seq)}, '
                f'asked for region {start+1}-{end})')
            return None
        left_padding = min(padding, start)
        right_padding = min(padding, len(seq) - end)
        return left_padding, right_padding, seq[start - left_padding : end + right_padding]


def load_sequences(filename, sequence, tag):
    if filename is not None and sequence is not None:
        panic(f'Both {tag} file and sequence were provided')
    if filename is None and sequence is None:
        panic(f'Neither {tag} file nor sequence were provided')

    res = Sequences()
    if sequence is not None:
        res.sequence = sequence.upper()
        return res

    if os.path.exists(filename + '.fai'):
        try:
            fasta = pysam.FastaFile(filename)
            res.indexed_fasta = fasta
        except (ValueError, OSError) as e:
            pass

    if res.indexed_fasta is None:
        res.sequences = {}
        with pysam.FastxFile(filename) as fastx:
            for record in fastx:
                res.sequences[record.name] = record.sequence.upper()
    return res


OPERATION_REGEX = re.compile(r'([0-9]+)([MIDNSHP=X])')

def process_cigar(split_cigar):
    """
    Processes CIGAR and returns a list of pairs (query_pos, ref_pos),
    as well as the alignment query and reference lengths.
    """
    ref_pos = 0
    query_pos = 0
    pairs = []
    for length, operation in split_cigar:
        length = int(length)
        match operation:
            case 'M' | '=' | 'X':
                for i in range(length):
                    pairs.append((query_pos + i, ref_pos + i))
                ref_pos += length
                query_pos += length
            case 'I' | 'S':
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
            case 'H':
                logging.error('Hard clipping may not be fully supported')
    return pairs, query_pos, ref_pos


LETTER_COLORS = { 'A': Fore.green, 'C': Fore.blue, 'G': Fore.yellow, 'T': Fore.red, '-': Fore.light_gray }
OP_COLORS = { ' ': '', '|': '', 'X': Fore.red, '-': Fore.red }


def print_cigar(split_cigar, out, painter, margin, active_width):
    s = painter.new_string()
    shift = 0
    for length, op in split_cigar:
        if shift == 0:
            s += margin
        elif shift + len(length) + len(op) >= active_width:
            shift = 0
            s += '\n' + margin
        else:
            shift += 1
            s += ' '
        s.extend(length, '')
        s.extend(op, Fore.green if op == 'M' or op == '=' else Fore.red)
        shift += len(length) + len(op)
    s += '\n'
    out.write(str(s))


def visualize_alignment(aln, ref_seqs, query_seqs, out, painter, padding, width, highlight, short_header):
    r_padd_left, r_padd_right, ref_seq = ref_seqs.fetch(aln.ref_name, aln.ref_start, aln.ref_end, padding)
    q_padd_left, q_padd_right, query_seq = (0, 0, aln.query_seq) if aln.query_seq is not None \
        else query_seqs.fetch(aln.query_name, aln.query_start, aln.query_end, padding)
    if ref_seq is None or query_seq is None:
        out.write(painter(f'\n\n!!! Unavailable sequences "{aln.query_name}" and "{aln.ref_name}" !!!\n\n', Fore.red))
        return
    if not aln.query_strand:
        query_seq = reverse_complement(query_seq)
        q_padd_left, q_padd_right = q_padd_right, q_padd_left

    split_cigar = OPERATION_REGEX.findall(aln.cigar)
    padded_cigar = deque(split_cigar)
    if q_padd_left:
        padded_cigar.appendleft((str(q_padd_left), 'I'))
    if r_padd_left:
        padded_cigar.appendleft((str(r_padd_left), 'D'))
    if r_padd_right:
        padded_cigar.append((str(r_padd_right), 'D'))
    if q_padd_right:
        padded_cigar.append((str(q_padd_right), 'I'))

    pairs, query_len, ref_len = process_cigar(padded_cigar)
    if len(ref_seq) != ref_len or len(query_seq) != query_len:
        sum_r_padd = r_padd_left + r_padd_right
        sum_q_padd = q_padd_left + q_padd_right
        logging.error(f'Alignment {aln.query_name} to {aln.ref_name}: incorrect alignment length: '
            f'reference {ref_len - sum_r_padd} (expected {len(ref_seq) - sum_r_padd}), '
            f'query {query_len - sum_q_padd} (expected {len(query_seq) - sum_q_padd})')
        out.write(painter(f'\n\n!!! Incorrect alignment {aln.query_name} to {aln.ref_name} !!!\n\n', Fore.red))
        return

    out.write('=====\n')
    out.write(painter(f'{aln.query_name}:{aln.query_start+1:,}-{aln.query_start + query_len:,}', Style.bold))
    out.write(' (top) to ')
    out.write(painter(f'{aln.ref_name}:{aln.ref_start+1:,}-{aln.ref_start + ref_len:,}', Style.bold))
    out.write(' (bottom)\n')

    if not short_header:
        print_cigar(split_cigar, out, painter, '  ', width - 2)

    hl_equal = 'E' in highlight
    hl_mism = 'X' in highlight
    hl_indels = 'I' in highlight

    ref_start = aln.ref_start - r_padd_left
    query_start = (aln.query_start if aln.query_strand else -(aln.query_end - 1)) - q_padd_left

    for i in range(0, len(pairs), width):
        qpos_line = ''
        rpos_line = ''
        qdigit_line = painter.new_string()
        rdigit_line = painter.new_string()
        aln_line = painter.new_string()
        qletter_line = painter.new_string()
        rletter_line = painter.new_string()

        for j, (query_offset, ref_offset) in enumerate(pairs[i:i + width]):
            if query_offset is None:
                qdigit_line += ' '
                qletter = '-'
            else:
                qpos = abs(query_start + query_offset) + 1
                digit = qpos % 10
                qdigit_line.extend(digit, '' if digit == 0 else Fore.light_gray)
                qletter = query_seq[query_offset]
                if (qpos == 1 or qpos % 50 == 0) and len(qpos_line) + 1 < j:
                    qpos_comma = f'{qpos:,}'
                    if j + 1 + len(qpos_comma) < width:
                        qpos_line += f'{" " * (j - len(qpos_line))}┌{qpos_comma}'

            if ref_offset is None:
                rdigit_line += ' '
                rletter = '-'
            else:
                rpos = ref_start + ref_offset + 1
                digit = rpos % 10
                rdigit_line.extend(digit, '' if digit == 0 else Fore.light_gray)
                rletter = ref_seq[ref_offset]
                if (rpos == 1 or rpos % 50 == 0) and len(rpos_line) + 1 < j:
                    rpos_comma = f'{rpos:,}'
                    if j + 1 + len(rpos_comma) < width:
                        rpos_line += f'{" " * (j - len(rpos_line))}└{rpos_comma}'

            if query_offset is None or ref_offset is None:
                op = '-' if hl_indels else ' '
            elif qletter != rletter:
                op = 'X' if hl_mism else ' '
            else:
                op = '|' if hl_equal else ' '

            aln_line.extend(op, OP_COLORS[op])
            qletter_line.extend(qletter, LETTER_COLORS.get(qletter, ''))
            rletter_line.extend(rletter, LETTER_COLORS.get(rletter, ''))

        out.write('\n')
        out.write(f'{qpos_line}\n')
        out.write(f'{qdigit_line}\n')
        out.write(f'{qletter_line}\n')
        out.write(f'{aln_line}\n')
        out.write(f'{rletter_line}\n')
        out.write(f'{rdigit_line}\n')
        out.write(f'{rpos_line}\n')


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
        help='Input SAM/BAM/CRAM file. Mutually exclusive with -p and -c.\n'
            'Does not require query sequences (-Q/-q).')
    aln_group.add_argument('-p', '--paf', metavar='FILE',
        help='Input PAF file. Mutually exclusive with -s and -c.')
    aln_group.add_argument('-c', '--cigar', metavar='STR',
        help='Cigar string. Mutually exclusive with -s and -p.\n'
            'Region definitions can be provided with -r/-q.')

    seq_args = parser.add_argument_group('Sequences')
    seq_args.add_argument('-R', '--ref-fa', metavar='FILE',
        help='Optionally indexed FASTA/Q file with reference sequences.')
    seq_args.add_argument('-Q', '--query-fa', metavar='FILE',
        help='Optionally indexed FASTA/Q file with query sequences.\n'
            'Defaults to -R if -Q/-q are not provided.')
    seq_args.add_argument('-r', '--ref-seq', metavar='STR',
        help='Reference sequence or coordinates (name:start-end, requires -R).')
    seq_args.add_argument('-q', '--query-seq', metavar='STR',
        help='Query sequence or coordinates (name:start-end[:+/-], requires -R or -Q).')

    opt_args = parser.add_argument_group('Optional arguments')
    opt_args.add_argument('-o', '--output', metavar='FILE', default='-',
        help='Output file [stdout]. Unless piped, stdout is automatically paginated.')
    opt_args.add_argument('-P', '--padding', metavar='INT', type=int, default=0,
        help='Add INT bp padding around the alignment region [%(default)s].')
    opt_args.add_argument('-H', '--highlight', metavar='STR', default='EX',
        help='Highlight these operations between two sequences [%(default)s], a string with the\n'
            'following possible letters: E (equal nucleotides), X (mismatches), I (indels).')
    opt_args.add_argument('-w', '--width', metavar='INT', default='auto',
        help='Display width [%(default)s].')
    opt_args.add_argument('-C', '--color', metavar='STR', default='auto',
        help='When to provide colored output (always|yes|y / auto / never|n) [%(default)s].')
    opt_args.add_argument('-S', '--short-header', action='store_true',
        help='Do not print CIGAR before alignment.')
    opt_args.add_argument('-h', '--help', action='help',
        help='Show this help message and exit.')
    opt_args.add_argument('-V', '--version', action='version',
        help="Show program's version number and exit.",
        version='visualize-cigar {version}'.format(version=importlib.metadata.version('visualize-cigar')))
    args = parser.parse_args()

    logger = logging.basicConfig(format='%(levelname)-7s| %(message)s', level=logging.DEBUG)

    force_ref_region = None
    force_query_region = None
    force_query_strand = None
    if args.ref_seq is not None and ':' in args.ref_seq:
        if args.ref_fa is None:
            panic('Fasta -R required since -r contains coordinates')
        if args.cigar is None:
            panic('Coordinates in -r are supported only together with -c CIGAR')
        force_ref_region = Interval.parse(args.ref_seq)
        args.ref_seq = None
    if args.query_seq is not None and ':' in args.query_seq:
        if args.ref_fa is None and args.query_fa is None:
            panic('Fasta -R/-Q required since -q contains coordinates')
        if args.cigar is None:
            panic('Coordinates in -q are supported only together with -c CIGAR')
        if args.query_seq.endswith(':-') or args.query_seq.endswith(':+'):
            force_query_strand = args.query_seq[-1] == '+'
            args.query_seq = args.query_seq[:-2]
        force_query_region = Interval.parse(args.query_seq)
        args.query_seq = None

    ref_seqs = load_sequences(args.ref_fa, args.ref_seq, 'reference')
    if args.sam is not None:
        query_seqs = None
    elif args.query_fa is not None or args.query_seq is not None:
        query_seqs = load_sequences(args.query_fa, args.query_seq, 'query')
    elif args.ref_fa is not None:
        query_seqs = ref_seqs
    else:
        panic('At least one of -R/-Q/-q arguments must be provided.')

    args.padding = max(0, args.padding)

    with open(args.output, 'wt', True) as out:
        if args.width == 'auto':
            args.width = os.get_terminal_size()[0] - 5 if out.isatty() else 120
        else:
            args.width = int(args.width)
        if args.width < 20:
            args.width = 20
            logging.warning(f'Width too small, expanding to {args.width}')

        painter = Painter(args.color, out)
        for aln in load_alignments(args, force_ref_region, force_query_region, force_query_strand):
            visualize_alignment(aln, ref_seqs, query_seqs, out, painter,
                args.padding, args.width, args.highlight, args.short_header)


if __name__ == '__main__':
    try:
        main()
    except (BrokenPipeError, KeyboardInterrupt) as e:
        pass
