import gzip
from PROFIT_seq.utils import to_str, to_bytes


def revcomp(seq):
    """
    Convert sequence to reverse complementary
    """
    trantab = str.maketrans("ATCG", "TAGC")
    return seq.translate(trantab)[::-1]


def is_gz_file(fname):
    """
    Detect magic bytes for gzip files
    From: https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed
    """
    with open(fname, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


def is_fq_file(fname):
    """
    Detect fastq/fasta
    """
    f = gzip.open(fname, 'rb') if is_gz_file(fname) else open(fname, 'r')
    header = to_str(f.readline())[0]
    f.close()
    return header == '@'


def yield_fq(fh):
    for line in fh:
        read_id = to_str(line).rstrip().lstrip('@')
        seq = to_str(fh.readline()).rstrip()
        sep = to_str(fh.readline()).rstrip()
        qual = to_str(fh.readline()).rstrip()
        yield (read_id, seq, sep, qual)


def yield_fa(fh):
    seq_id = None
    seq = ''
    for line in fh:
        tmp_line = to_str(line)
        if tmp_line.startswith('>'):
            if seq_id is not None:
                yield (seq_id, seq, None, None)
            seq_id = tmp_line.rstrip().lstrip('>')
            seq = ''
        else:
            seq += tmp_line.rstrip()
    yield (seq_id, seq, None, None)


def yield_fastx(fname):
    """
    Wrapper for iter return fasta record
    """
    is_gz = is_gz_file(fname)
    is_fq = is_fq_file(fname)

    fh = gzip.open(fname, 'rb') if is_gz else open(fname, 'r')
    if is_fq:
        yield from yield_fq(fh)
    else:
        yield from yield_fa(fh)


def load_fastx(fname):
    """
    Wrapper for load fasta/fastq record
    """
    sequences = {}
    for seq_id, seq, sep, qual in yield_fastx(fname):
        sequences[seq_id] = seq
    return sequences


def count_fastx(fname):
    return sum([1 for _ in yield_fastx(fname)])
