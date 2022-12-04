from y_dc3 import dc3
import cProfile


def bwt(string, end_of_string="$"):
    """
    Compute the BWT from the suffix

    Args:
        string (str): string
        end_of_string (char): appended character specifying the end of
            the string

    Return:
        bwt (str): BWT transformation of T
    """
    bwtStr = ""
    string += end_of_string

    s_table = dc3(string, is_seq_to_ascii=True)
    print(f"Suffix table (dc3): {s_table}")

    for suffix_pos in s_table:
        bwtStr += string[suffix_pos - 1]

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


def create_rank_mat(dnaSeq: str) -> dict:
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
    lenStr = len(dnaSeq)
    for letter in alphabet:
        # :(i + 1) because we want to include the first & last char of string
        rkMat[letter] = [dnaSeq[:(i + 1)].count(letter) for i in range(lenStr)]
    return rkMat


def search_kmer_pos(genome: str, kmer: str):
    """Search a pattern in a String using the BWT.

    Args:
        genome (str): a string made only with A, T, C, G
        suffix_tab (list[int]): sorted array of all the suffixes of
            genome
        pattern (str): pattern we are looking for

    Return:
        bool: true if the pattern is in the string
    """
    isKmerIn = False

    # To be passed in argument to avoid computing each time
    bwtGen = bwt(genome)
    bwtSort = list(bwtGen)
    bwtSort.sort()

    lenBwt = len(bwtGen)
    e = 0
    f = lenBwt - 1  # Stay in the string boundaries
    i = len(kmer) - 1  # Stay in the kmer boundaries

    rankMat = create_rank_mat(bwtGen)
    print(f"Rank matrix: {rankMat}")
    nbOccur = 1

    while nbOccur > 0 and i >= 0:
        X = kmer[i]

        if X not in genome:
            return False
        else:
            print(f"i: {i}, X: {X}")
            print(f"Rank table of {X}: {rankMat[X]}")

            firstOcc = bwtSort.index(X) - 1  # $ isn't taken into account
            print(f"First index of {X} in bwtSort: {firstOcc}")

            # The first character of the BWT has a rank of 1 in the rank
            # matrix but it is the first appearance of this character.
            # To take into account this information, we decrease e of 1 in
            # presence of this character, when there is no character before it
            if len(bwtGen[:(e + firstOcc)]) == 0:
                e = firstOcc + rankMat[X][e] - 1
            else:
                e = firstOcc + rankMat[X][e]
            f = firstOcc + rankMat[X][f]
            print(f"e: {e}, f: {f}")

            nbOccur = f - e  # quantity of elements between 2 indexes
            if nbOccur == 0:
                return False
            print(f"Number of pattern: {nbOccur}\n")

            i -= 1

            # True if the entire pattern is crossed
            if i == -1:
                isKmerIn = True
                return isKmerIn, nbOccur

    return isKmerIn, nbOccur


if __name__ == "__main__":
    T = "ATAATAGCCGAAAGCCTTTAGTACGCCAATGGTTAACACACA" * 300000
    cProfile.run("dc3(T,)")

    # print(f"Rank matrix of T: {create_rank_mat(bwtT)}")  # OK

    # print(search_kmer_pos(T, "ATA"))
