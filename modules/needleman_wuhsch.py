import pandas
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from Bio.Align import substitution_matrices
from itertools import product


""" 
Cost functions for the Needleman-Wunsch algorithm: of 2 symbols, of n+1 symbols, of n+m symbols
"""

def cost_2_symbols_alignment(i: chr, j: chr, blosum_m: bool, identity_score=1, substitution_score=-1):
    """
    Calculate the score of aligning 2 characters: i and j. If blosum_m is True, we use BLOSUM62 matrix.

    Parameters:
    ----------
    i : chr
        First character to align.
    j : chr
        Second character to align.
    blosum_m : bool
        If True, we use BLOSUM62 matrix.
    identity_score : int
        Score for aligning identical characters.
    substitution_score : int
        Score for aligning non-identical characters.

    Returns:
    -------
    score : int
        The score of aligning i and j.
    """
    if blosum_m:
        matrix = substitution_matrices.load("BLOSUM62")
        return int(matrix[(i, j)])

    if i == j:
        return identity_score
    else:
        return substitution_score

def cost_n_and_1_alignment(i, alignements, blosum_m: bool, gap_score, identity_score=1, substitution_score=-1):
    """
    Calculate the score of aligning character i with list of characters written in alignements. If blosum_m is True, we use BLOSUM62 matrix.

    Parameters:
    ----------
    i : chr
        One character to align with others.
    alignements : list of chr
        Already aligned characters.
    blosum_m : bool
        If True, we use BLOSUM62 matrix.
    gap_score : int
        The score of opening a gap. 
    identity_score : int
        Score for aligning identical characters.
    substitution_score : int
        Score for aligning non-identical characters.

    Returns:
    -------
    score : int
        The score of aligning i and j.
    """
    score = 0
    if i == '-':
        return gap_score
    else:
        for c in alignements:
            if c == '-':
                score += gap_score
            else:
                score += cost_2_symbols_alignment(i, c, blosum_m, identity_score, substitution_score)
    return score / len(alignements)

def cost_n_and_m_alignment(list1, list2, blosum_m=False, identity_score=1, substitution_score=-1):
    """ 
    returns the cost of aligning two character lists
    Character lists represent the already aligned characters of the blocks
    Parameters:
    ------------
    list1 : list of already aligned caracters 1
    list2 : list of already aligned caracters 2
    blosum_m : bool
        If True, we use BLOSUM62 matrix.
    identity_score : int
        Score for aligning identical characters.
    substitution_score : int
        Score for aligning non-identical characters.

    Returns:
    Score : int, associated score
    """
    score=0
    n=len(list1)
    m=len(list2)

    if blosum_m:
        matrix = substitution_matrices.load("BLOSUM62")

    for i in range(n):
        for j in range(m):
            if blosum_m:
                if list1[i]!='-' and list2[j]!='-':
                    score+= int(matrix[(list1[i],list2[j])])
                if list1[i]=='-' and list2[j]=='-':
                    score+=identity_score
                else :
                    score+=0
            else :
                if list1[i] == list2[j]:
                    score+= identity_score
                else:
                    score+= substitution_score
    return score/(n*m)


"""
Plotting and printing functions for the Needleman-Wunsch algorithm
"""

def plot_nw_matrix(matrix, arrow_matrix, block1, block2):
    """
    Visualize the Needleman-Wunsch matrix with arrows for n+m sequences.

    Parameters:
    ----------
    matrix : pandas.DataFrame
        The filled matrix with scores.
    arrow_matrix : pandas.DataFrame
        The matrix of arrows indicating traceback paths.
    block1 : list of str
        List containing sequences already aligned
    block 2 : list of str
        List containing sequences already aligned. Objective is now to align block1 and block2
    """
    fig, ax = plt.subplots(figsize=(6, 4))
    
    # Convert matrix to numpy array, replacing None with 0
    matrix_np = matrix.fillna(0).to_numpy()
    
    # Define labels for the axes (insert '-' at the beginning to account for initial gap)
    # row_labels = ['-'] + [ [block1[k][i] for k in range (len(block1))] for i in range (len(block1[0])) ]
    # col_labels = ['-'] + [ [block2[k][i] for k in range (len(block2))] for i in range (len(block2[0])) ]

    # Define labels for the axes (insert '-' at the beginning to account for initial gap)
    if isinstance(block1, str):
        block1 = [block1]
    if isinstance(block2, str):
        block2 = [block2]

    row_labels = ['-'] + list(block1[0])
    col_labels = ['-'] + list(block2[0])

    for i in range(1, len(col_labels)):
        for j in range(1, min(len(block2), 6)):
            col_labels[i] = col_labels[i] + '\n' + block2[j][i-1]
    
    for i in range(1, len(row_labels)):
        for j in range(1, min(len(block1), 6)):
            row_labels[i] = row_labels[i] + '\n' + block1[j][i-1]
    
    # Create a heatmap with custom labels
    
    sns.heatmap(matrix_np.astype(float), annot=True, fmt=".0f", cmap="Blues", linewidths=0.5, 
                ax=ax, cbar=False, xticklabels=col_labels, yticklabels=row_labels)

    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")

    traceback_indexes = [(matrix.shape[0]-1, matrix.shape[1]-1)]
    while traceback_indexes[-1] != (0, 0):
        for prev_i, prev_j in arrow_matrix.at[traceback_indexes[-1][0], traceback_indexes[-1][1]]:
            traceback_indexes.append((prev_i, prev_j))

    # Draw arrows based on arrow_matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if arrow_matrix.at[i, j] is not None:
                for prev_i, prev_j in arrow_matrix.at[i, j]:
                    dx = prev_j - j  # X direction
                    dy = prev_i - i  # Y direction

                    # color is red if arrow is in traceback_indexes
                    if ((prev_i, prev_j) in traceback_indexes and (i, j) in traceback_indexes):
                        color = 'red'
                    else:
                        color = 'black'
                    if dx == 0:
                        ax.arrow(j + 0.5, i + 0.35, dx * 0.5, dy * 0.5, head_width=0.1, 
                                 head_length=0.1, fc=color, ec=color)
                    elif dy == 0:
                        ax.arrow(j + 0.25, i + 0.5, dx * 0.4, dy * 0.5, head_width=0.1, 
                                 head_length=0.1, fc=color, ec=color)
                    else:
                        ax.arrow(j + 0.35, i + 0.35, dx * 0.6, dy * 0.6, head_width=0.1, 
                             head_length=0.1, fc=color, ec=color)
    
    ax.set_title("Needleman-Wunsch Alignment Matrix with Traceback Arrows")
    plt.show()

def plot_nw_matrix_3_seq(matrix, sequences):
    '''
    Plot the Needleman-Wunsch matrix for 3 sequences.

    Parameters:
    ----------
    matrix : dict
        Dictionary where keys are index tuples (i1, i2, i3) and values are DP scores.
    sequences : list of str
        List of 3 sequences to align.
    '''
    matrix_dict = matrix

    sequences = [f"-{seq}" for seq in sequences]

    # Dimensions of the matrix
    max_i = max(k[0] for k in matrix_dict) + 1
    max_j = max(k[1] for k in matrix_dict) + 1
    max_k = max(k[2] for k in matrix_dict) + 1


    matrix_np = np.full((max_i, max_j, max_k), np.nan)  
    for (i, j, k), value in matrix_dict.items():
        matrix_np[i, j, k] = value

    # Plotting the matrix
    fig, axes = plt.subplots(1, max_k, figsize=(5 * max_k, 5))

    if max_k == 1:
        axes = [axes]  

    for k in range(max_k):
        ax = axes[k]
        sns.heatmap(matrix_np[:, :, k], annot=True, cmap="coolwarm", fmt=".1f", linewidths=0.5, ax=ax,
                    xticklabels=list(sequences[1]), yticklabels=list(sequences[0])) 
        
        ax.xaxis.set_ticks_position("top") 
        ax.xaxis.set_label_position("top")  
        ax.set_title(f"Slice at {sequences[2][k]}", fontsize=12) 
  
    fig.suptitle("Needleman-Wunsch Matrix N^4", fontsize=16, fontweight="bold")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()


def print_alignments(alignments):
    """
    Print the aligned sequences.

    Parameters:
    ----------
    alignments : list of str
        List containing the aligned sequences.
    """
    print("Alignments:")
    for alignment in alignments:
        print(alignment)
    return

    
def print_nw_result_multidim_3_seq(matrix, sequences, score, alignments):
    '''
    Print the Needleman-Wunsch result for multidimensional case with 3 sequences.

    Parameters:
    ----------
    matrix : pandas.DataFrame
        The filled matrix.
    sequences : list of str
        List of 3 sequences to align.
    score : int
        The score of the alignment.
    alignments : list of str
        The aligned sequences.
    '''
    print(f"Alignment was made with Multidimensional Needleman-Wunsch algorithm. Score is {score}.")
    print_alignments(alignments)
    plot_nw_matrix_3_seq(matrix, sequences)
    return


def print_nw_result_multidim_N_seq(score, alignments):
    '''
    Print the Needleman-Wunsch result for multidimensional case with 3 sequences.

    Parameters:
    ----------
    matrix : pandas.DataFrame
        The filled matrix.
    sequences : list of str
        List of 3 sequences to align.
    score : int
        The score of the alignment.
    alignments : list of str
        The aligned sequences.
    '''
    print(f"Alignment was made with Multidimensional Needleman-Wunsch algorithm. Score is {score}.")
    print_alignments(alignments)
    return


def print_nw_result(matrix, arrow_matrix, score, alignments, block1, block2):
    '''
    Print the Needleman-Wunsch result.

    Parameters:
    ----------
    matrix : pandas.DataFrame
        The filled matrix.
    arrow_matrix : pandas.DataFrame
        The matrix of arrows. In matrix_arrows[i, j] we store the indexes of the cells where we can go from cell i, j.
    score : int
        The score of the alignment.
    alignments : list of str
        The aligned sequences.
    block1 : list of str
        List containing sequences already aligned
    block 2 : list of str
        List containing sequences already aligned. Objective is now to align block1 and block2
    '''
    print(f"Alignment was made with Needleman-Wunsch algorithm. Score is {score}.")
    print_alignments(alignments)
    plot_nw_matrix(matrix, arrow_matrix, block1, block2)
    return


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
    Fill the Needleman-Wunsch matrix up to N^K and store the indexes of the arrows with possibility to have up to K arrows.

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
        Dictionary where keys are index tuples (i1, ..., iK) and values store previous indices for traceback.
    '''

    K = len(sequences)

    # Determine the matrix dimensions for each axis (N_i + 1)
    shape = [len(seq) + 1 for seq in sequences]

    # Initialize DP-matrix, traceback and gap-matrices
    matrix = {}
    arrow_matrix = {}
    gap_matrix = {}

    # Step 1: Initialize DP matrix (boundaries)
    for index in product(*[range(n) for n in shape]):
        if sum(index) == 0:
            matrix[index] = 0
            arrow_matrix[index] = None
            gap_matrix[index] = 0
        else:
            min_gap = min(i for i in index if i > 0)
            matrix[index] = gap_opening_score + (min_gap - 1) * gap_extension_score
            arrow_matrix[index] = []
            gap_matrix[index] = 1

    # Step 2: Filling the DP matrix
    for index in product(*[range(n) for n in shape]):
        if sum(index) == 0:
            continue  # (0,0,...,0)

        best_score = float('-inf')
        best_prev = None
        best_gap_state = None

        # We go through all possible steps (inserting gaps in different sequences)
        for shift in product([0, -1], repeat=K):
            if sum(shift) == 0:
                continue  # (0,0,...,0)

            prev_index = tuple(i + s for i, s in zip(index, shift))
            if any(i < 0 for i in prev_index):
                continue  # Ignore negative indices

            # Form a list of characters to be aligned
            aligned_chars = [
                sequences[i][index[i] - 1] if index[i] > 0 else '-'
                for i in range(K)
            ]

            gap_state = gap_matrix[prev_index] == 1  # Check if we are in a gap state
            gap_penalty = gap_extension_score if gap_state else gap_opening_score

            # Calculate the cost of aligning the last character
            cost = cost_n_and_1_alignment(
                aligned_chars[-1],  
                aligned_chars[:-1],  
                blosum_m,
                gap_penalty,
                identity_score,
                substitution_score
            )

            score = matrix[prev_index] + cost

            # Update the best score and previous index
            if score > best_score:
                best_score = score
                best_prev = prev_index
                best_gap_state = 1 if '-' in aligned_chars else 0 

        # Store the best score and previous index
        matrix[index] = best_score
        arrow_matrix[index] = best_prev
        gap_matrix[index] = best_gap_state

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