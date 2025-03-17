from Bio.Align import substitution_matrices

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

def cost_n_symbols_alignment(symbols, blosum_m: bool, gap_penalty, identity_score=1, substitution_score=-1):
    """
    Computes the cost of aligning a set of symbols (of dimension K) in the multidimensional Niedelman-Wunsch algorithm.

    Parameters:
    ----------
    symbols : list of str
        List of characters that are aligned at the same time (например, ['A', 'G', 'A'] или ['C', '-', 'T']).
    blosum_m : bool
        If True, then the BLOSUM62 matrix is used to calculate the score.
    gap_penalty : int
    identity_score : int
    substitution_score : int

    Returns:
    -------
    score : int
        Cost of aligning the symbols.
    """
    if '-' in symbols:
        return gap_penalty

    if blosum_m:
        matrix = substitution_matrices.load("BLOSUM62")
        score = 0
        for i in range(len(symbols)):
            for j in range(i + 1, len(symbols)):  # Take only unique pairs
                score += int(matrix[(symbols[i], symbols[j])])
        return score / (len(symbols) * (len(symbols) - 1) / 2)  # Mean score of pairs

    score = 0
    for i in range(len(symbols)):
        for j in range(i + 1, len(symbols)):
            if symbols[i] == symbols[j]:
                score += identity_score
            else:
                score += substitution_score

    return score / (len(symbols) * (len(symbols) - 1) / 2) 