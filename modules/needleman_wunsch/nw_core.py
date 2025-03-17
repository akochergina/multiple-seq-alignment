import pandas
import pandas as pd
from itertools import product
import sys
import os

# Add the parent directory to the path so we can import the module
sys.path.append(os.path.abspath(os.path.join('..')))

from modules.needleman_wunsch.nw_cost import *
from modules.needleman_wunsch.nw_visualization import *
from modules.carillo_lipman_bounds import *

"""
Needleman-Wunsch functions.
"""

def fill_needleman_wunsch_matrix(sequence, previous_alignment, blosum_m, gap_opening_score, gap_extension_score, identity_score=1, substitution_score=-1):
    '''
    Fill the Needleman-Wunsch matrix and store the indexes of the arrows with possibility to have up to 3 arrows.

    Parameters:
    ----------
    sequence :  str
        Sequences to align with the previous alignment.
    previous_alignment : list of str
        Some sequences already aligned.
    blosum_m : bool
        If True, we use BLOSUM62 matrix.
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
    matrix : pandas.DataFrame
        The filled matrix.
    arrow_matrix : pandas.DataFrame
        The matrix of arrows indicating traceback paths.
    '''
    # Initialize the matrix
    rows = len(sequence) + 1
    cols = len(previous_alignment[0]) + 1

    matrix = pandas.DataFrame(index=range(rows), columns=range(cols))
    arrow_matrix = pandas.DataFrame(index=range(rows), columns=range(cols))
    gap_matrix = pandas.DataFrame(index=range(rows), columns=range(cols))

    arrow_matrix.at[0, 0] = None
    matrix.at[0, 0] = 0
    gap_matrix.at[0, 0] = 0
 

    # Initialize the first row and column
    for i in range(1, rows):
        matrix.at[i, 0] = gap_opening_score + (i-1) * gap_extension_score
        arrow_matrix.at[i, 0] = [(i-1, 0)]
        gap_matrix.at[i, 0] = 1
    for j in range(1, cols):
        matrix.at[0, j] = gap_opening_score + (j-1) * gap_extension_score
        arrow_matrix.at[0, j] = [(0, j-1)]
        gap_matrix.at[0, j] = 1


    # Fill the matrix 
    for i in range(1, rows):
        for j in range(1, cols):
            # Calculate the scores
            chr_alignment = [previous_alignment[k][j-1] for k in range(len(previous_alignment))]
            if blosum_m:
                if gap_matrix.at[i-1, j] == 1:
                    match = matrix.at[i-1, j-1] + cost_n_and_1_alignment(sequence[i-1], chr_alignment, blosum_m, gap_extension_score)
                else:
                    match = matrix.at[i-1, j-1] + cost_n_and_1_alignment(sequence[i-1], chr_alignment, blosum_m, gap_opening_score)
            else:
                if gap_matrix.at[i-1, j] == 1:
                    match = matrix.at[i-1, j-1] + cost_n_and_1_alignment(sequence[i-1], chr_alignment, blosum_m, gap_extension_score, identity_score, substitution_score)
                else:
                    match = matrix.at[i-1, j-1] + cost_n_and_1_alignment(sequence[i-1], chr_alignment, blosum_m, gap_opening_score, identity_score, substitution_score)
            if gap_matrix.at[i-1, j] == 1:
                delete = matrix.at[i-1, j] + gap_extension_score
            else:
                delete = matrix.at[i-1, j] + gap_opening_score
            if gap_matrix.at[i, j-1] == 1:
                insert = matrix.at[i, j-1] + gap_extension_score
            else:
                insert = matrix.at[i, j-1] + gap_opening_score

            # Update the matrix
            max_score = max(match, delete, insert)
            arrow_matrix.at[i, j] = []
            matrix.at[i, j] = max_score
            if max_score == match:
                arrow_matrix.at[i, j].append((i-1, j-1))
                gap_matrix.at[i, j] = 0
            if max_score == delete:
                arrow_matrix.at[i, j].append((i-1, j))
                gap_matrix.at[i, j] = 1
            if max_score == insert:
                arrow_matrix.at[i, j].append((i, j-1))
                gap_matrix.at[i, j] = 1

    return matrix, arrow_matrix

def fill_needleman_wunsch_matrix_multiple(block1, block2, blosum_m, gap_opening_score, gap_extension_score, identity_score=1, substitution_score=-1):
    '''
    Fill the Needleman-Wunsch matrix for multiple sequences, regrouped in two already aligned blocks n+m
    and store the indexes of the arrows with possibility to have up to 3 arrows.

    Parameters:
    ----------
    block1 : list of str
        List containing sequences already aligned
    block 2 : list of str
        List containing sequences already aligned. Objective is now to align block1 and block2
    blosum_m : bool
        If True, we use BLOSUM62 matrix.
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
    matrix : pandas.DataFrame
        The filled matrix.
    '''
    # Initialize the matrix
    rows = len(block1[0]) + 1
    cols = len(block2[0]) + 1
    matrix = pandas.DataFrame(index=range(rows), columns=range(cols))
    arrow_matrix = pandas.DataFrame(index=range(rows), columns=range(cols))
    gap_matrix = pandas.DataFrame(index=range(rows), columns=range(cols))

    arrow_matrix.at[0, 0] = None
    matrix.at[0, 0] = 0
    gap_matrix.at[0, 0] = 0
 

    # Initialize the first row and column
    for i in range(1, rows):
        matrix.at[i, 0] = gap_opening_score + (i-1) * gap_extension_score
        arrow_matrix.at[i, 0] = [(i-1, 0)]
        gap_matrix.at[i, 0] = 1
    for j in range(1, cols):
        matrix.at[0, j] = gap_opening_score + (j-1) * gap_extension_score
        arrow_matrix.at[0, j] = [(0, j-1)]
        gap_matrix.at[0, j] = 1

    # Fill the matrix 
    for i in range(1, rows):
        for j in range(1, cols):
            # Calculate the scores
            match = matrix.at[i-1, j-1] + cost_n_and_m_alignment([(block1[k][i-1]) for k in range(len(block1))],[(block2[k][j-1]) for k in range (len(block2))], blosum_m, identity_score, substitution_score)
            if gap_matrix.at[i-1, j] == 1:
                delete = matrix.at[i-1, j] + gap_extension_score
            else:
                delete = matrix.at[i-1, j] + gap_opening_score
            if gap_matrix.at[i, j-1] == 1:
                insert = matrix.at[i, j-1] + gap_extension_score
            else:
                insert = matrix.at[i, j-1] + gap_opening_score

            # Update the matrix
            max_score = max(match, delete, insert)
            arrow_matrix.at[i, j] = []
            matrix.at[i, j] = max_score
            if max_score == match:
                arrow_matrix.at[i, j].append((i-1, j-1))
                gap_matrix.at[i, j] = 0
            if max_score == delete:
                arrow_matrix.at[i, j].append((i-1, j))
                gap_matrix.at[i, j] = 1
            if max_score == insert:
                arrow_matrix.at[i, j].append((i, j-1))
                gap_matrix.at[i, j] = 1

    return matrix, arrow_matrix

def fill_needleman_wunsch_matrix_multidim(sequences, blosum_m, gap_opening_score, gap_extension_score, identity_score=1, substitution_score=-1):
    '''
    Fill the Needleman-Wunsch matrix up to N^K and store the index of the best previous step (not multiple paths).

    Parameters:
    ----------
    sequences : list of str
        List of K sequences to align.
    blosum_m : bool
        If True, we use BLOSUM62 matrix.
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
    matrix : dict
        Dictionary where keys are index tuples (i1, ..., iK) and values are DP scores.
    arrow_matrix : dict
        Dictionary where keys are index tuples (i1, ..., iK) and values store previous index for traceback.
    '''
    K = len(sequences)
    shape = [len(seq) + 1 for seq in sequences]

    matrix = {}
    arrow_matrix = {}
    gap_matrix = {}

    # Step 1: Initialize DP matrix within bounds
    for index in product(*[range(n) for n in shape]):
        matrix[index] = 0 if sum(index) == 0 else sum(gap_opening_score + (i - 1) * gap_extension_score for i in index if i > 0)
        arrow_matrix[index] = None
        gap_matrix[index] = 0 if sum(index) == 0 else 1
    
    # Step 2: Fill DP matrix within bounds
    for index in product(*[range(n) for n in shape]):
        if sum(index) == 0:
            continue  

        score_dict = {}  # Store scores for each possible previous index
        best_prev, best_score, best_gap_state = None, float('-inf'), None
        
        for shift in product([0, -1], repeat=K):
            if sum(shift) == 0:
                continue  

            prev_index = tuple(i + s for i, s in zip(index, shift))
            if any(i < 0 for i in prev_index):
                continue   # Skip if out of bounds

            aligned_chars = []
            for i in range(K):
                if prev_index[i] == index[i]:  # Если индекс не изменился, значит, тут разрыв
                    aligned_chars.append('-')
                else:
                    aligned_chars.append(sequences[i][index[i] - 1])

            gap_state = gap_matrix[prev_index] == 1
            gap_penalty = gap_extension_score if gap_state else gap_opening_score

            cost = cost_n_symbols_alignment(
                aligned_chars,  
                blosum_m,
                gap_penalty,
                identity_score,
                substitution_score
            )
        

            score = matrix[prev_index] + cost
            score_dict[prev_index] = score  # Store computed score


            if score > best_score:
                best_prev, best_score = prev_index, score
                best_gap_state = 1 if '-' in aligned_chars else 0
    

        matrix[index] = best_score
        arrow_matrix[index] = best_prev
        gap_matrix[index] = best_gap_state
        
        # print(f"Index: {index}, Scores: {score_dict}, Best: {best_prev} -> {best_score}")

    return matrix, arrow_matrix


def fill_needleman_wunsch_matrix_multidim_carillo(sequences, blosum_m, gap_opening_score, gap_extension_score, identity_score=1, substitution_score=-1):
    """
    Fill the Needleman-Wunsch matrix using Carillo-Lipman bounds and measure efficiency.

    Parameters:
    ----------
    sequences : list of str
        List of K sequences to align.
    blosum_m : bool
        If True, we use BLOSUM62 matrix.
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
    matrix : dict
        Dictionary where keys are index tuples (i1, ..., iK) and values are DP scores.
    arrow_matrix : dict
        Dictionary where keys are index tuples (i1, ..., iK) and values store previous index for traceback.
    """

    K = len(sequences)
    shape = [len(seq) + 1 for seq in sequences]

    # Compute Carillo-Lipman bounds
    bounds = compute_carillo_lipman_bounds(sequences, blosum_m, True, gap_opening_score, gap_extension_score, identity_score, substitution_score)

    # Compute max possible gaps for each sequence (before the last letter!)
    max_gap_seq = tuple(bounds[i][1] - len(sequences[i]) for i in range(K))

    valid_moves = False

    matrix = {}
    arrow_matrix = {}
    gap_matrix = {}
    available_gaps = {}  # Dict to track remaining gaps for each sequence

    # Анализ Carillo-Lipman: считаем заполненные и пропущенные ячейки
    total_possible_cells = np.prod(shape)  # Максимально возможное количество ячеек
    filled_cells = 0  # Заполненные ячейки
    skipped_cells = 0  # Пропущенные ячейки (Carillo-Lipman)

    # Step 1: Initialize DP matrix within bounds
    for index in product(*[range(n) for n in shape]):
        if any(index[i] < bounds[i][0] or index[i] > bounds[i][1] for i in range(K)):
            skipped_cells += 1
            continue 

        filled_cells += 1
        matrix[index] = 0 if sum(index) == 0 else sum(gap_opening_score + (i - 1) * gap_extension_score for i in index if i > 0)
        arrow_matrix[index] = None
        gap_matrix[index] = 0 if sum(index) == 0 else 1
        available_gaps[index] = max_gap_seq  # Initialize available gaps for each index

    # Step 2: Fill DP matrix within bounds
    for index in product(*[range(n) for n in shape]):
        if sum(index) == 0:
            continue

        if any(index[i] < bounds[i][0] or index[i] > bounds[i][1] for i in range(K)):
            continue  

        score_dict = {}  # Store scores for each possible previous index
        best_prev, best_score, best_gap_state = None, float('-inf'), None
        best_gaps = None  # Track best available gaps

        for shift in product([0, -1], repeat=K):
            if sum(shift) == 0:
                continue  

            prev_index = tuple(i + s for i, s in zip(index, shift))

            if any(prev_index[i] < bounds[i][0] or prev_index[i] > bounds[i][1] for i in range(K)):
                continue

            if any(i < 0 for i in prev_index):
                continue  
            
            # Проверяем, не станет ли `available_gaps` отрицательным
            new_gaps = list(available_gaps.get(prev_index, max_gap_seq))
            for i in range(K):
                if prev_index[i] == index[i] and index[i] < len(sequences[i]):  
                    new_gaps[i] -= 1

            if any(gap < 0 for gap in new_gaps):
                continue  # Если хотя бы один гэп < 0, этот путь невозможен

            valid_moves = True  # Нашли хотя бы один допустимый `prev_index`
            new_gaps = tuple(new_gaps)

            aligned_chars = []
            for i in range(K):
                if prev_index[i] == index[i]:  
                    aligned_chars.append('-')
                else:
                    aligned_chars.append(sequences[i][index[i] - 1])

            gap_state = gap_matrix.get(prev_index, 0) == 1
            gap_penalty = gap_extension_score if gap_state else gap_opening_score

            cost = cost_n_symbols_alignment(
                aligned_chars,  
                blosum_m,
                gap_penalty,
                identity_score,
                substitution_score
            )

            score = matrix.get(prev_index, float('-inf')) + cost
            score_dict[prev_index] = score 

            # Update available gaps
            new_gaps = list(available_gaps.get(prev_index, max_gap_seq))
            for i in range(K):
                if prev_index[i] == index[i] and index[i] < len(sequences[i]):  
                    new_gaps[i] -= 1
            new_gaps = tuple(new_gaps) 

            if score > best_score:
                best_prev, best_score = prev_index, score
                best_gap_state = 1 if '-' in aligned_chars else 0
                best_gaps = new_gaps
        
        if not valid_moves:
            matrix[index] = float('-inf')
            arrow_matrix[index] = None
            available_gaps[index] = max_gap_seq
            skipped_cells += 1
        else:
            matrix[index] = best_score
            arrow_matrix[index] = best_prev
            gap_matrix[index] = best_gap_state
            available_gaps[index] = best_gaps if best_gaps else max_gap_seq  # Store updated gap info



    # Results Carillo-Lipman
    print()
    efficiency = (skipped_cells) / total_possible_cells if total_possible_cells > 0 else 0
    print(f"Total possible cells: {total_possible_cells}")
    print(f"Skipped cells due to Carillo-Lipman bounds: {skipped_cells}")
    print(f"Efficiency of Carillo-Lipman bounds: {efficiency:.2%}")

    return matrix, arrow_matrix


def needleman_wunsch_step(sequence, previous_alignment, blosum_m, gap_opening_score=-10, gap_extension_score=-2, print_result=False, identity_score=1, substitution_score=-1):
    """
    Perform Needleman-Wunsch alignment of 1 sequence with the previous alignment of N.

    Parameters:
    ----------
    sequence : str
        One sequence to align with the previous alignment.
    previous_alignment : list of str
        Some sequences already aligned.
    blosum_m : bool
        If True, we use BLOSUM62 matrix.
    gap_opening_score : int
        Score for opening a gap.
    gap_extension_score : int
        Score for extending a gap.
    print_result : bool
        If True, print the matrix.
    identity_score : int
        Score for aligning identical characters.
    substitution_score : int
        Score for aligning non-identical characters.

    Returns:
    -------
    alignment : tuple
        A tuple containing the aligned sequences and a score.
    """
    
    if blosum_m:
        matrix, arrow_matrix = fill_needleman_wunsch_matrix(sequence, previous_alignment, blosum_m, gap_opening_score, gap_extension_score)
    else:
        matrix, arrow_matrix = fill_needleman_wunsch_matrix(sequence, previous_alignment, blosum_m, gap_opening_score, gap_extension_score, identity_score, substitution_score)
    
    score = matrix.at[len(sequence), len(previous_alignment[0])]


    alignement1 = ''
    alignement2 = ''
    i = len(sequence)
    j = len(previous_alignment[0])

    while i > 0 or j > 0:
        prev_i, prev_j = arrow_matrix.at[i, j][0]
        if i - prev_i == 1 and j - prev_j == 1:
            alignement1 = sequence[i-1] + alignement1
            alignement2 = previous_alignment[0][j-1] + alignement2
        elif i - prev_i == 1:
            alignement1 = sequence[i-1] + alignement1
            alignement2 = '-' + alignement2
        else:
            alignement1 = '-' + alignement1
            alignement2 = previous_alignment[0][j-1] + alignement2
        i = prev_i
        j = prev_j

    sequences = [sequence] + previous_alignment

    if len(alignement1) != len(previous_alignment[0]):
        max_length = max(len(alignement1), len(alignement2))
        if len(alignement1) < max_length:
            alignement1 = alignement1.ljust(max_length, '-')
        else:
            previous_alignment = [seq.ljust(max_length, '-') for seq in previous_alignment]


    if len(previous_alignment) < 2:
        previous_alignment = [alignement2]

    previous_alignment.insert(0, alignement1)

    if print_result:
        print_nw_result(matrix, arrow_matrix, score, previous_alignment, sequences[0], sequences[1:])

    return score, previous_alignment

def needleman_wunsch(sequences, blosum_m, gap_opening_score=-10, gap_extension_score=-2, print_result=False, identity_score=1, substitution_score=-1):
    """
    Perform Needleman-Wunsch alignment.

    Parameters:
    ----------
    sequences : list of str
        Sequences to align.
    blosum_m : bool
        If True, we use BLOSUM62 matrix.
    gap_opening_score : int
        Score for opening a gap.
    gap_extension_score : int
        Score for extending a gap.
    print_result : bool
        If True, print the matrix.
    identity_score : int
        Score for aligning identical characters.
    substitution_score : int
        Score for aligning non-identical characters.

    Returns:
    -------
    alignment : tuple
        A tuple containing the aligned sequences and a score.
    """
    previous_score, previous_alignment = needleman_wunsch_step(sequences[1], [sequences[0]], blosum_m, gap_opening_score, gap_extension_score, print_result, identity_score, substitution_score)


    for i in range(2, len(sequences)):
        previous_score, previous_alignment = needleman_wunsch_step(sequences[i], previous_alignment, blosum_m, gap_opening_score, gap_extension_score, print_result, identity_score, substitution_score)
    
    return previous_score, previous_alignment


def needleman_wunsch_multiple(block1, block2, blosum_m, gap_opening_score=-10, gap_extension_score=-2, print_result=False, identity_score=1, substitution_score=-1):
    """
    Perform Needleman-Wunsch alignment, to align two blocks of already aligned sequences

    Parameters:
    ----------
    block1 : list of str
        List containing sequences already aligned
    block 2 : list of str
        List containing sequences already aligned. Objective is now to align block1 and block2
    blosum_m : bool
        If True, we use BLOSUM62 matrix.
    gap_opening_score : int
        Score for opening a gap.
    gap_extension_score : int
        Score for extending a gap.
    print_result : bool
        If True, print the matrix.
    identity_score : int
        Score for aligning identical characters.
    substitution_score : int
        Score for aligning non-identical characters.

    Returns:
    -------
    alignment : tuple
        A tuple containing the aligned sequences and a score.
    """
    i = len(block1[0])
    j = len(block2[0])
    if i < j :
        block1, block2 = block2, block1
        i, j = j, i

    if blosum_m:
        matrix, arrow_matrix = fill_needleman_wunsch_matrix_multiple(block1, block2, blosum_m, gap_opening_score, gap_extension_score)
    else:
        matrix, arrow_matrix = fill_needleman_wunsch_matrix_multiple(block1, block2, blosum_m, gap_opening_score, gap_extension_score, identity_score, substitution_score)
    
    score = matrix.at[len(block1[0]), len(block2[0])]


    aa = [("") for i in range (len(block1))]
    bb = [("") for i in range (len(block2))]

    while i > 0 or j > 0:
        prev_i, prev_j = arrow_matrix.at[i, j][0]
        if i - prev_i == 1 and j - prev_j == 1:
            for k in range(len(block1)):
                if i > 0:  
                    aa[k] += block1[k][i-1]  
                else:
                    aa[k] += '-'
            for k in range(len(block2)):
                if j > 0:
                    bb[k] += block2[k][j-1]
                else:
                    bb[k] += '-'
        elif i - prev_i == 1:
            for k in range (len(block1)):
                aa[k]+=(block1[k][i-1])
            for k in range(len(block2)):
                bb[k]+='-'
        else:
            for k in range (len(block1)):
                aa[k]+='-'
            for k in range(len(block2)):
                bb[k]+=(block2[k][j-1])
        i = prev_i
        j = prev_j
    
    
    for k in range (len(block1)):
      aa[k]=aa[k][::-1]
    for k in range(len(block2)):
      bb[k]=bb[k][::-1]
    alignments = aa+bb

    if print_result:
        print_nw_result(matrix, arrow_matrix, score, alignments, block1, block2)

    return score, alignments


def needleman_wunsch_multidim(sequences, blosum_m, carillo, gap_opening_score=-10, gap_extension_score=-2, print_result=False, identity_score=1, substitution_score=-1):
    """
    Perform Needleman-Wunsch alignment for multiple sequences using a multidimensional DP matrix.

    Parameters:
    ----------
    sequences : list of str
        List of K sequences to align.
    blosum_m : bool
        If True, use BLOSUM62 matrix.
    carillo : bool
        If True, use Carillo-Lipman bounds.
    gap_opening_score : int
        Score for opening a gap.
    gap_extension_score : int
        Score for extending a gap.
    print_result : bool
        If True, print the result.
    identity_score : int
        Score for aligning identical characters.
    substitution_score : int
        Score for aligning non-identical characters.

    Returns:
    -------
    score : int
        The optimal alignment score.
    alignments : list of str
        The aligned sequences.
    """
    if len(sequences) == 2:
        return needleman_wunsch(sequences, blosum_m, gap_opening_score, gap_extension_score, print_result, identity_score, substitution_score)

    # Step 1: Fill the DP matrix
    if carillo:
        matrix, arrow_matrix = fill_needleman_wunsch_matrix_multidim_carillo(
            sequences, blosum_m, gap_opening_score, gap_extension_score, identity_score, substitution_score
        )
    else:
        matrix, arrow_matrix = fill_needleman_wunsch_matrix_multidim(
            sequences, blosum_m, gap_opening_score, gap_extension_score, identity_score, substitution_score
        )

    # Step 2: Backtracking to reconstruct the alignment
    K = len(sequences)
    shape = [len(seq) for seq in sequences]
    index = tuple(n for n in shape)  # Start at the bottom-right corner


    score = matrix[index]
    aligned_sequences = [""] * K

    while index is not None and sum(index) > 0:
        prev_index = arrow_matrix[index]  # Only one best path

        # Build aligned sequences
        for i in range(K):
            if prev_index and prev_index[i] == index[i] - 1:
                aligned_sequences[i] = sequences[i][index[i] - 1] + aligned_sequences[i]
            else:
                aligned_sequences[i] = "-" + aligned_sequences[i]

        index = prev_index  # Move to the previous index

    # Step 3: Print the results if needed
    if print_result:
        if K == 2:
            print_nw_result(matrix, arrow_matrix, score, aligned_sequences, sequences[0], sequences[1:])
        elif K == 3:
            print_nw_result_multidim_3_seq(matrix, sequences, score, aligned_sequences)
        else:
            print_nw_result_multidim_N_seq(score, aligned_sequences)

    return score, aligned_sequences
