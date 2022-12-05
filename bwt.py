import numpy as np
from copy import deepcopy


def bwt(dna, suffixTable, end_of_string="$"):
    """Returns the Burrows Wheeler Transform of the dna argument.
    The BWT is computed efficiently using the suffix array of dna.

    Args:
        dna (str): string made only with A, T, C and G character
        suffix_table (ndarray): contains all the suffixes of dna
        end_of_string (char): appended character specifying the end of
            the string

    Return:
        str: BWT transformation of dna
    """
    bwtStr = ""
    dna += end_of_string

    for suffix_pos in suffixTable:
        bwtStr += dna[suffix_pos - 1]

    return bwtStr


def create_rank_mat(bwtDna):
    """Returns a dict storing the rank tables of each character in
    the dnaSeq argument.
    A rank table is an array containing, for each character of bwtDna
    and at each position, the number of occurences of the character
    until now.

    Args:
        dnaSeq (str): a string made only with A, T, C, G and a $

    Returns:
        dict: contains 5 keys: A, T, C, G and $, each key is mapped
            with the rank table in the dnaSeq of the associated
            character
    """
    alphabet = ['A', 'T', 'C', 'G', '$']
    rkMat = {}
    lenDna = len(bwtDna)
    for letter in alphabet:
        # We use uint32 because 32 bits is sufficient to store the
        # unsigned values we need
        countings = np.zeros(lenDna, dtype=np.uint32)
        # We count only one over 32 rows to be time efficient
        for i in range(lenDna):
            if i % 32 == 0:
                # :(i + 1) because we want to include the first & last char of string
                countings[i] = bwtDna[:(i + 1)].count(letter)
            else:
                countings[i] = 0
        rkMat[letter] = countings
    return rkMat


def search_kmer_pos(bwtDna, kmer, suffixTab):
    """Search a kmer in a Burrows Wheeler tranformed genome. Returns
    the kmer and it's position(s) in the genome.

    Args:
        bwtDna (str): the BWT of a string made only with A, T, C, G
            and a $
        kmer (str): nucleotides pattern to search in the genome
        suffixTab (ndarray): sorted array of all the suffixes
            used to built bwtDna

    Return:
        bool: true if the pattern is in the string
    """
    isKmerIn = False

    bwtSort = deepcopy(bwtDna).sort()  # Avoid modifying bwtGen while sorting
    lenBwt = len(bwtDna)

    e = 0
    f = lenBwt - 1  # Stay in the string boundaries
    i = len(kmer) - 1  # Stay in the kmer boundaries

    rankMat = create_rank_mat(bwtDna)
    # print(f"Rank matrix: {rankMat}")
    nbOccur = 1

    while nbOccur > 0 and i >= 0:
        X = kmer[i]

        if X not in bwtDna:
            return False
        else:
            # print(f"i: {i}, X: {X}")
            # print(f"Rank table of {X}: {rankMat[X]}")

            firstOcc = bwtSort.index(X) - 1  # $ isn't taken into account
            # print(f"First index of {X} in bwtSort: {firstOcc}")

            newE = e
            rank = rankMat[X][newE]
            counts = 0
            while rank == 0:
                newE -= 1
                if bwtDna[newE] == X:
                    counts -= 1
                rank = rankMat[X][newE]
            rank -= counts

            # TODO Do the same but for f

            # The first character of the BWT has a rank of 1 in the
            # rank matrix but it is the first appearance of this
            # character. To take into account this information, we
            # decrease e of 1 in presence of this character, i.e when
            # there is no character before it
            # TODO integrate this part in the while above
            if len(bwtDna[:(e + firstOcc)]) == 0:
                e = firstOcc + rankMat[X][e] - 1
            else:
                e = firstOcc + rankMat[X][e]
            f = firstOcc + rankMat[X][f]
            # print(f"e: {e}, f: {f}")

            nbOccur = f - e  # quantity of elements between 2 indexes
            if nbOccur == 0:
                return False
            # print(f"Number of pattern: {nbOccur}\n")

            i -= 1

            # True if the entire pattern is crossed
            if i == -1:
                isKmerIn = True
                return isKmerIn, nbOccur

    return isKmerIn, nbOccur


if __name__ == "__main__":
    T = "ATAATAGGATCCGA" * 500
    rankMat = create_rank_mat(T)
    disp = "[ "
    for count in rankMat['A']:
        disp += f" {count}"
    disp += "]"
    print(disp)

    # bwtT = bwt(T)
    # print(f"BWT from suffix table: {bwtT}")

    # bwtDC3 = bwt(T)
    # print(f"Is the BWT DC3 the same ? {bwtDC3 == bwtT}")
    # print(f"BWT DC3: {bwtDC3}")

    # print(search_kmer_pos(T, "ATA"))
