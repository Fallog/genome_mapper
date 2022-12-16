import numpy as np


def bwt(dna: str, suffix_t: np.ndarray, string_end: str = "$") -> str:
    """Returns the Burrows Wheeler Transform of the dna argument.
    The BWT is computed efficiently using the suffix array of dna.

    Args:
        dna (str): string made only with A, T, C and G character
        suffix_t (np.ndarray): contains all the suffixes of dna, output
            of dc3 algorithm over dna
        string_end (str): appended character specifying the end of
            the string

    Return:
        str: Burrows Wheeler Transform of dna
    """
    bwt = ""

    for i in range(len(suffix_t)):
        # Chromosomes have sometimes bases in lowercase
        if suffix_t[i] == 0:
            bwt += string_end
        else:
            bwt += dna[suffix_t[i] - 1].capitalize()

    return bwt


def create_rank_mat(bwt_dna: str) -> dict[np.ndarray]:
    """Returns a dict storing the rank tables of each character in
    the dnaSeq argument.
    A rank table is an array containing, for each character of bwt_dna
    and at each position, the number of occurences of the character
    until this position.

    Args:
        bwt_dna (str): a string made only with A, T, C, G and a $

    Returns:
        dict[np.ndarray]: contains 5 keys: A, T, C, G and $, each key is
        mapped with its rank table
    """
    alphabet = ["A", "T", "C", "G", "$"]
    dna_len = len(bwt_dna)
    rank_mat = {
        "A": np.zeros(dna_len + 1, dtype=np.int32),
        "T": np.zeros(dna_len + 1, dtype=np.int32),
        "C": np.zeros(dna_len + 1, dtype=np.int32),
        "G": np.zeros(dna_len + 1, dtype=np.int32),
        "$": np.zeros(dna_len + 1, dtype=np.int32),
    }
    for i, char in enumerate("o" + bwt_dna):
        if i == 0:
            for other_letter in alphabet:
                rank_mat[other_letter][i] = 0
        else:
            for other_letter in alphabet:
                rank_mat[other_letter][i] = rank_mat[other_letter][i - 1]
            rank_mat[char][i] += 1
    return rank_mat
