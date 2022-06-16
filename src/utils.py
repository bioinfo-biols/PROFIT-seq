import ssw
from seqIO import revcomp


def grouper(iterable, n, fillvalue=None):
    """
    Collect data info fixed-length chunks or blocks
    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    """
    from itertools import zip_longest
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=None)


def align_sequence(aligner, query, strandness=True):
    for_aln = aligner.align(query)
    if strandness:
        return for_aln

    rev_aln = aligner.align(revcomp(query))
    if for_aln.score > rev_aln.score:
        return for_aln
    else:
        return rev_aln



