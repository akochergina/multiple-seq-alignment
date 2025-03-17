import sys
import os
from Bio import pairwise2
from Bio.Align import substitution_matrices
from itertools import permutations

# Add the parent directory to the path so we can import the module
sys.path.append(os.path.abspath(os.path.join('..')))

from modules.needleman_wunsch.nw_core import *
from modules.needleman_wunsch.nw_visualization import *



def sp_score(alignment1, alignment2):
    """
    Computes Sum-of-Pairs (SP) score to compare two multiple sequence alignments.

    Parameters:
    ----------
    alignment1 : list of str
        The first alignment (e.g., from user's algorithm).
    alignment2 : list of str
        The reference alignment (e.g., from BALIbase).

    Returns:
    -------
    float
        The similarity score (1 = identical, 0 = completely different).
    """
    total_pairs = 0
    matching_pairs = 0

    for col in zip(*alignment1):
        total_pairs += 1
        if col in zip(*alignment2):  # Check if column exists in reference alignment
            matching_pairs += 1

    return matching_pairs / total_pairs if total_pairs > 0 else 0


def alignment_result_comparison(result1, result2):
    """ 
    Compare two alignment results, considering sequence order variations.
    
    Parameters:
    ----------
    result1 : tuple (float, list of str, list of str)
        Custom Needleman-Wunsch alignment result in the form (score, alignment1, alignment2).
    result2 : tuple (float, list of str, list of str)
        Biopython Needleman-Wunsch alignment result in the same format.
    
    Prints:
    ------
    - Whether the scores match.
    - The highest similarity between alignments using the sum-of-pairs (SP) score.
    """

    score1, alignment1_1, alignment1_2 = result1
    score2, alignment2_1, alignment2_2 = result2

    print(f"Custom Needleman-Wunsch Score: {score1}")
    print(f"Biopython Needleman-Wunsch Score: {score2}")

    if score1 == score2:
        print("✅ The alignment scores match!")
    else:
        print("❌ The alignment scores differ.")

    # Compute SP-score for both orderings and take the maximum
    sp_similarity = max(
        sp_score([alignment1_1, alignment1_2], [alignment2_1, alignment2_2]),
        sp_score([alignment1_1, alignment1_2], [alignment2_2, alignment2_1])
    )

    print(f"SP Similarity Score: {sp_similarity:.4f}")

    if sp_similarity == 1.0:
        print("✅ Alignments are identical.")
    elif sp_similarity > 0.8:
        print("⚠️ Alignments are highly similar but not identical.")
    else:
        print("❌ Alignments differ significantly.")


def test_biopython_vs_custom(seq1, seq2, blossum=True, gap_open_penalty=-10, gap_extension_penalty=-2):
    substitution_matrix = substitution_matrices.load("BLOSUM62") if blossum else substitution_matrices.load("DNAFULL")

    alignments_biopython = pairwise2.align.globalds(seq1, seq2, substitution_matrix, gap_open_penalty, gap_extension_penalty)
    print("Biopython results:")
    print_alignments([alignments_biopython[0][0], alignments_biopython[0][1]])
    print(" \nCustom results:")
    score, alignment = needleman_wunsch([seq1, seq2], True, print_result=False)
    print_alignments(alignment)
    print()

    custom_result = (score, alignment[0], alignment[1])
    result_biopython = (alignments_biopython[0][2], alignments_biopython[0][0], alignments_biopython[0][1])
    alignment_result_comparison(custom_result, result_biopython)


def alignment_result_comparison_3seq(result1, result2):
    """ 
    Compare two alignment results for three sequences, considering all sequence order variations.
    
    Parameters:
    ----------
    result1 : tuple (float, list of str)
        Custom Needleman-Wunsch alignment result in the form (score, alignment1, alignment2, alignment3).
    result2 : tuple (float, list of str)
        Another Needleman-Wunsch alignment result in the same format.
    
    Prints:
    ------
    - Whether the scores match.
    - The highest similarity between alignments using the sum-of-pairs (SP) score.
    """

    score1, alignment1_1, alignment1_2, alignment1_3 = result1
    score2, alignment2_1, alignment2_2, alignment2_3 = result2

    print(f"Custom Needleman-Wunsch Score: {score1}")
    print(f"Biopython Needleman-Wunsch Score: {score2}")

    if score1 == score2:
        print("✅ The alignment scores match!")
    else:
        print("❌ The alignment scores differ.")

    # Generate all permutations of alignments for comparison
    custom_alignments = [alignment1_1, alignment1_2, alignment1_3]
    other_alignments = [alignment2_1, alignment2_2, alignment2_3]

    max_sp_similarity = max(
        sp_score(custom_alignments, list(permuted_alignments))
        for permuted_alignments in permutations(other_alignments)
    )

    print(f"SP Similarity Score: {max_sp_similarity:.4f}")

    if max_sp_similarity == 1.0:
        print("✅ Alignments are identical.")
    elif max_sp_similarity > 0.8:
        print("⚠️ Alignments are highly similar but not identical.")
    else:
        print("❌ Alignments differ significantly.")