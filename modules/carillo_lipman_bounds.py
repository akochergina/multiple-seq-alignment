import numpy as np
from itertools import combinations
import sys
import os

# Add the parent directory to the path so we can import the module
sys.path.append(os.path.abspath(os.path.join('..')))

from modules.needleman_wunsch.nw_core import *
from modules.needleman_wunsch.nw_visualization import *

def seq_min_max(sequence):
    """
    Compute the minimum and maximum index of a sequence, ignoring leading and trailing '-'.
    
    Parameters:
    ----------
    sequence : str
        The aligned sequence containing characters and '-'.
    
    Returns:
    -------
    seq_min : int
        Index of the first non '-' character.
    seq_max : int
        Index of the last non '-' character.
    """
    seq_min = next((i for i, char in enumerate(sequence) if char != '-'), None)
    seq_max = next((i for i, char in enumerate(reversed(sequence)) if char != '-'), None)
    
    if seq_min is None:
        return None, None  # Sequence contains only '-'
    
    seq_max = len(sequence) - 1 - seq_max
    
    return seq_min, seq_max

def print_carillo_lipman_bounds(bounds, sequences):
    """
    Print the Carillo-Lipman bounds.

    Parameters:
    ----------
    bounds : list of tuples
        List of tuples (min_i, max_i) for each sequence.
    sequences : list of str
        List of K sequences to align.
    """
    K = len(sequences)
    print("Carillo-Lipman bounds:")
    for i in range(K):
        print(f"Sequence {i}: {bounds[i]} - {sequences[i]}")
    return

def compute_carillo_lipman_bounds(sequences, blosum_m, print_results=True, gap_opening_score=-10, gap_extension_score=-2, identity_score=1, substitution_score=-1):
    """
    Compute Carillo-Lipman bounds by performing pairwise alignments.

    Parameters:
    ----------
    sequences : list of str
        List of K sequences to align.
    blosum_m : bool
        If True, use BLOSUM62 matrix.
    gap_opening_score : int
        Score for opening a gap.
    gap_extension_score : int
        Score for extending a gap.
    identity_score : int
        Score for aligning identical characters.
    substitution_score : int
        Score for aligning non-identical characters.

    Returns:
    -------
    bounds : list of tuples
        List of tuples (min_i, max_i) for each sequence.
    """

    K = len(sequences)
    pairwise_bounds = {i: [] for i in range(K)}

    # 1. Compute pairwise alignments
    for i, j in combinations(range(K), 2):
        _, aligned = needleman_wunsch([sequences[j], sequences[i]], blosum_m, gap_opening_score, gap_extension_score, identity_score=identity_score, substitution_score=substitution_score)

        # Find min and max index mapping for each sequence
        seq_i = aligned[0]
        seq_j = aligned[1]

        min_i, max_i = seq_min_max(seq_i)
        min_j, max_j = seq_min_max(seq_j)

        pairwise_bounds[i].append((min_i, max_i))
        pairwise_bounds[j].append((min_j, max_j))

        print_alignments(aligned)
        print(f"Sequence {i}: {min_i} - {max_i}")
        print(f"Sequence {j}: {min_j} - {max_j}")
        print()

            

    # 2. Compute global bounds for each sequence
    bounds = []
    for i in range(K):
        min_bound = min(b[0] for b in pairwise_bounds[i])
        max_bound = max(b[1] for b in pairwise_bounds[i])
        bounds.append((min_bound, max_bound))
    
    bounds = [(min_bound, max_bound + 1) for min_bound, max_bound in bounds]

    if print_results:
        print_carillo_lipman_bounds(bounds, sequences)

    return bounds