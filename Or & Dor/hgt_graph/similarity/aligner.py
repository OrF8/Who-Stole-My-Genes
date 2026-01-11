from Bio.Align import PairwiseAligner
from ..constants import BLOSUM62, GAP_EXTEND_PENALTY, GAP_OPEN_PENALTY

def make_aligner() -> PairwiseAligner:
    """
    Create and configure a PairwiseAligner for local sequence alignment using BLOSUM62 matrix.
    :return: Configured PairwiseAligner object.
    """
    aligner: PairwiseAligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.substitution_matrix = BLOSUM62
    aligner.open_gap_score = GAP_OPEN_PENALTY
    aligner.extend_gap_score = GAP_EXTEND_PENALTY
    return aligner
