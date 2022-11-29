from y_dc3 import dc3, np


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


def bwt(T, end_of_string="$"):
    """
    Compute the BWT from the suffix table

    Args:
        T (str): string
        end_of_string (char): end of string character to append

    Return:rankXInSubL
        bwt (str): BWT transformation of T
    """
    bwtStr = ""

    T += end_of_string
    s_table = suffix_table(T)  # Has to be replace by DC3 algorithm

    for tuple in s_table:
        index = tuple[1]
        bwtStr += T[index - 1]
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


def create_rank_table(string: str) -> np.array:
    """Returns the array of the rank of each character in the string
    argument. The rank of a character represents the number of
    occurences of the same character before it in the string.

    Args:
        string (str): a string

    Returns:
        np.array[int]: contains the ranks of each character in string
    """
    strLen = len(string)
    rank = np.empty(strLen, dtype=int)  # Performance
    for i in range(strLen):
        # Counting elements before bwt[i] that are similar to bwt[i]
        rank[i] = string[:i].count(string[i])

    return rank


def search_kmer_pos(genome: str, kmer: str):
    """Search a pattern in a String using the BWT

    Args:
        S (str): string
        pattern (str): pattern we are looking for

    Return:
        bool: true if the pattern is in the string
    """
    isKmerIn = False

    bwtGen = bwt(genome)
    lenBwt = len(bwtGen)
    bwtSort = list(bwtGen)
    bwtSort.sort()
    e = 0
    f = lenBwt
    i = len(kmer) - 1
    rank = create_rank_table(bwtGen)
    nbOccur = 1

    while nbOccur > 0 and i >= 0:
        X = kmer[i]
        print(f"i: {i}, X: {X}")

        print(f"e: {e}, f: {f}")

        if X not in genome:
            return False
        else:
            # Because the slicing stops on the character before f, we add the
            # remaining character to not lose any information.
            subL = bwtGen[e:f]
            if f < lenBwt:
                subL += bwtGen[f]

            print(f"subL: {subL}")

            # If f is the size of F, we don't need to append the last character
            # because it is included in the slicing of subRank
            print(f"rank[e:f]: {rank[e:f]}")
            if f < lenBwt:
                print(f"rank[f]: {rank[f]}, type: {type(rank[f])}")
                subRank = np.concatenate((rank[e:f], [rank[f]]))
            else:
                subRank = np.array(rank[e:f])
            # all of the ranks of characters that are present in subL
            print(f"subRank: {subRank}")

            # Create an empty array of the size of subRank and taking only
            # the filled part
            rankX = np.empty(subRank.size, dtype=int)

            j = 0  # increment for subRank going throug all the array
            k = 0  # increment for rankX stopping after
            while j < subRank.size:
                if subL[j] == X:
                    rankX[k] = subRank[j]
                    k += 1
                j += 1
            rankX = np.take(rankX, list(range(k)))
            print(f"rankXInSubL: {rankX}")

            if rankX.size == 0:
                nbOccur = 0
                break

            nbOccur = rankX.size

            frstOcc = bwtSort.index(X)

            # contains all the occurences of X in F
            e = frstOcc + rankX[0]
            f = frstOcc + rankX[-1]
            i -= 1

            # True if the entire pattern is crossed
            if i == -1:
                isKmerIn = True
                return isKmerIn, nbOccur

    return isKmerIn, nbOccur


if __name__ == "__main__":
    T = "abaaba"
    # T = "abcabcacab"
    arr = np.array([1, 3, 5, 7, 9])
    emptArr = np.array(arr[1])
    print(f"Test: {emptArr}, type: {type(emptArr)}")

    bwtT = bwt(T)
    print(bwtT)

    print(f"Rank table of T: {create_rank_table(T)}")

    print(f"Inverse BWT result: {efficient_inverse_BWT(bwtT)}")  # OK

    print(f"Original suffix table {suffix_table(T)}")
    # print(f"DC3 suffix table {y_dc3.dc3(T)}")  # TODO: À FAIRE MARCHER LOL
    print(create_rank_table(T) == [0, 0, 1, 2, 1, 3])  # OK

    print(search_kmer_pos(T, "ababa"))
