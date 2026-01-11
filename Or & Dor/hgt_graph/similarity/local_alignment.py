from typing import Tuple
from Bio.Align import PairwiseAligner

def local_alignment_metrics(aligner: PairwiseAligner, seq1: str, seq2: str) -> Tuple[float, float, int]:
    """
    Perform local alignment between two sequences and calculate identity, coverage, and aligned length.
    :param aligner: A configured PairwiseAligner object.
    :param seq1: The first sequence.
    :param seq2: The second sequence.
    :return: A tuple containing identity, coverage, and aligned length.
    """
    alignments = aligner.align(seq1, seq2)

    aln = alignments[0]

    aligned_len = 0
    matches = 0

    # aln.aligned gives aligned blocks as coordinate ranges
    for (s1_start, s1_end), (s2_start, s2_end) in zip(aln.aligned[0], aln.aligned[1]):
        block_len = s1_end - s1_start
        aligned_len += block_len

        for i in range(block_len):
            if seq1[s1_start + i] == seq2[s2_start + i]:
                matches += 1

    if aligned_len == 0:
        return 0.0, 0.0, 0

    identity = matches / aligned_len
    coverage = aligned_len / min(len(seq1), len(seq2))

    return identity, coverage, aligned_len