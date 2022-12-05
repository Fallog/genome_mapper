from y_dc3 import dc3, np
from copy import deepcopy
import array


def suffix_list(T):
    """
    Compute the suffix list of T argument

    Args:
        T (str): string

    Return:
        list of strings: suffix list
    """
    suffix_list = []
    for i in range(len(T)):
        suffix_list.append(T[(-i - 1):])
    return suffix_list


def suffix_table(T):
    """
    Compute the suffix table

    Args:
        T (str): string

    Return:
        list of tuples (suffix,location): suffix table, location being
            the index of the first character of the suffix in the total
            string
    """
    suffix_table = []
    suffix_table = suffix_list(T)
    i = 0
    j = len(suffix_table) - 1
    while i < len(suffix_table):
        suffix_table[i] = (suffix_table[i], j)
        i += 1
        j -= 1

    suffix_table.sort()
    return suffix_table


def bwt(string, end_of_string="$"):
    """
    Compute the BWT from the suffix table

    Args:
        string (str): string
        end_of_string (char): appended character specifying the end of
            the string

    Return:
        bwt (str): BWT transformation of T
    """
    bwtStr = ""

    string += end_of_string
    s_table = dc3(string)  # Has to be replace by DC3 algorithm

    for tuple in s_table:
        index = tuple[1]
        bwtStr += string[index - 1]
    return bwtStr


def efficient_inverse_BWT(bwtStr: str, end_of_string: str = "$") -> str:
    """
    Returns the original string that were used to build the bwtStr
    argument.

    Args:
        bwtStr (str): BWT of a string
        last_character (char): specifies the end of the string used in
            the bwt

    Return:
        T (str): BWT^{-1} of bwt without the last end_of_string char
    """
    T = ""
    lenStr = len(bwtStr)  # Performance
    appar_order_table = [0] * lenStr  # Performance
    for i in range(lenStr):
        # Counting elements before bwt[i] that are similar to bwt[i]
        appar_order_table[i] = 1 + bwtStr[:i].count(bwtStr[i])
    X = bwtStr[0]
    k = 1
    T += end_of_string

    while X != end_of_string:
        T = X + T
        j = k

        for i in range(len(bwtStr)):
            if bwtStr[i] < X:
                # Localisation of X in the sorted BWT, obtained by
                # counting every characters that are less than X in the BWT
                j += 1

        X = bwtStr[j - 1]
        k = appar_order_table[j - 1]

    return T[:-1]


def create_rank_mat(bwtDna: str) -> dict:
    """Returns a dict storing the rank tables of each character in
    the dnaSeq argument. Because the DNA is exclusively made of A, T,
    C and G, the returned dict will contain 5 keys, the 4 bases and a $
    specifying the end of the sequence.

    Args:
        dnaSeq (str): a string made only with A, T, C, G and a $

    Returns:
        dict: contains 5 character keys, each key is mapped with the
            rank table in the dnaSeq of the associated character
    """
    alphabet = ['A', 'T', 'C', 'G', '$']
    rkMat = {}
    lenDna = len(bwtDna)
    for letter in alphabet:
        # :(i + 1) because we want to include the first & last char of string
        countings = [0] * lenDna
        for i in range(lenDna):
            countings[i] = bwtDna[:(i + 1)].count(letter)
        rkMat[letter] = array.array('H', countings)
    return rkMat


def search_kmer_pos(bwtDna: str, kmer: str, suffixTab: np.ndarray) -> tuple[str, list[int]]:
    """Search a kmer in a Burrows Wheeler tranformed genome. Returns
    the kmer and it's position(s) in the genome.

    Args:
        bwtDna (str): the BWT of a string made only with A, T, C, G
            and a $
        kmer (str): nucleotides pattern to search in the genome
        suffixTab (np.ndarray[int]): sorted array of all the suffixes
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

            # The first character of the BWT has a rank of 1 in the
            # rank matrix but it is the first appearance of this
            # character. To take into account this information, we
            # decrease e of 1 in presence of this character, i.e when
            # there is no character before it
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

    bwtT = bwt(T)
    print(bwtT)

    print(f"Rank matrix of T: {create_rank_mat(bwtT)}")  # OK

    print(f"Inverse BWT result: {efficient_inverse_BWT(bwtT)}")  # OK

    print(search_kmer_pos(T, "ATA"))
