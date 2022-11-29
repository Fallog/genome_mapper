import y_dc3


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
        suffix_list.append(T[(- i - 1):])
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


def efficient_inverse_BWT(bwt, end_of_string="$"):
    """
    Inverse the BWT

    Args:
        bwt (str): bwt of a string T
        last_character (char): which is the end of string character?

    Return:
        T (str): BWT^{-1} of bwt
    """
    T = ""
    # Initialisation #
    appar_order_table = []
    for i in range(len(bwt)):
        # Counting elements before bwt[i] that are similar to bwt[i]
        appar_order_table.append(1 + bwt[:i].count(bwt[i]))
    X = bwt[0]
    k = 1
    T += end_of_string

    while X != end_of_string:
        T = X + T
        j = k

        for i in range(len(bwt)):
            if bwt[i] < X:
                # Localisation of X in the sorted BWT, obtained by
                # counting every characters that are less than X in the BWT
                j += 1

        X = bwt[j - 1]
        k = appar_order_table[j - 1]

    return (T[:-1])


def create_rank_table(string: str) -> list[int]:
    strLen = len(string)
    rank = [0] * strLen
    for i in range(strLen):
        # Counting elements before bwt[i] that are similar to bwt[i]
        rank[i] = string[:i].count(string[i])

    return rank


def search_kmer_pos(genome: str, kmer: str):
    """
    Search a pattern in a String using the BWT

    Args:
        S (str): string
        pattern (str): pattern we are looking for

    Return:
        bool: true if the pattern is in the string
    """
    isKmerIn = False
    L = bwt(genome)
    F = list(L)
    F.sort()
    e = 0
    f = len(F)
    i = len(kmer) - 1
    rank = create_rank_table(L)
    nbOccur = 1

    while nbOccur > 0 and i >= 0:
        X = kmer[i]

        # Because the slicing stops on the character before f, we add the
        # remaining character to not lose any information.
        subL = L[e:f] + L[f if f < len(F) else f:]
        print(f"subL: {subL}")

        # If f is the size of F, we don't need to append the last character
        # because it is included in the slicing of subRank
        # all of the ranks of characters that are present in subL
        subRank = rank[e:f] + ([rank[f]] if f < len(F) else [])
        print(f"subRank: {subRank}")

        rankXInSubL = []
        for j in range(len(subRank)):
            if subL[j] == X:
                rankXInSubL.append(subRank[j])
        print(f"rankList: {rankXInSubL}")

        if len(rankXInSubL) == 0:
            nbOccur = 0
            break

        nbOccur = len(rankXInSubL)

        if X not in subL:
            break
        else:
            fstOcc = F.index(X)
            e = F.index(X)  # e is the index of first X in the sorted BWT
            f = len(F) - 1 - F[::-1].index(X)  # same for f but with last X

            # contains all the occurences of X in F
            allX = F[e:f] + ([F[f]] if f < len(F) else [])
            print(f"All the occurences of {X}: {allX}")
            e = fstOcc + rankXInSubL[0]
            print(f"Index of the first occurence of {X} in F: {e}")
            f = fstOcc + rankXInSubL[-1]
            print(f"Index of the last occurence of {X} in F: {f}")

            i -= 1

    # False if the first letter is not in the pattern argument
    # True if the entire pattern is crossed
    if i == -1 and kmer[0] in genome:
        isKmerIn = True
        return isKmerIn, nbOccur

    return isKmerIn, nbOccur


if __name__ == "__main__":
    T = "abaaba"

    bwtT = bwt(T)
    print(bwtT)

print(f"Original suffix table {suffix_table(T)}")
# print(f"DC3 suffix table {y_dc3.dc3(T)}")  # TODO: À FAIRE MARCHER LOL
print(suffix_table(T))

print(create_rank_table(T) == [0, 0, 1, 2, 1, 3])  # OK

print(search_kmer_pos(T, "rudwaaaba"))
