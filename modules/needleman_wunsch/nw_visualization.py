import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

"""
Plotting and printing functions for the Needleman-Wunsch algorithm
"""

"""
1. Plotting matrix for different cases of Needleman-Wunsch algorithm
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


"""
2. Printing results for different cases of Needleman-Wunsch algorithm
"""

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
