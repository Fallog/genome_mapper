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

    Return:
        bwt (str): BWT transformation of T
    """
    bwt = ""

    T += end_of_string
    s_table = suffix_table(T)  # Has to be replace by DC3 algorithm

    for tuple in s_table:
        index = tuple[1]
        bwt += T[index - 1]
    return (bwt)


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


def pattern_matching_BWT(S, pattern):
    """
    Search a pattern in a String using the BWT

    Args:
        S (str): string
        pattern (str): pattern we are looking for

    Return:
        bool: true if the pattern is in the string
    """
    pattern_in_S = False
    L = bwt(S)
    F = list(L)
    F.sort()
    e = 1
    f = len(F)
    i = len(pattern) - 1

    while e < f and i > 0:
        X = pattern[i]

        # Extending the range of sub_L to the next character to avoid missing
        # a character and stop the search
        sub_L = L[e:f] + L[f if f < len(F) else len(F) - 1]

        if X not in sub_L:
            break
        else:
            e = F.index(X)  # e is the index of first X in the sorted BWT
            f = len(F) - 1 - F[::-1].index(X)  # same for f but with last X

            i -= 1

    # False if the first letter is not in pattern argument
    # True if the entire pattern is crossed
    if i == 0 and pattern[0] in S:
        pattern_in_S = True
    return pattern_in_S


if __name__ == "__main__":
    T = "ACATACAGATG"

    bwtT = bwt(T)
    print(bwtT)

    print(suffix_table(T))
