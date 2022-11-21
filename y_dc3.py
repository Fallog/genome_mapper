# import numpy   ## A convertir en numpy
from array import array
import numpy as np
from quickSort import quickSort


def strToAscii(text: str) -> array:  # ya un numpy loadtxt qui existe
    """Transform a string to is Ascii value ant put it in a numpy array
    TO DO : update for DNA (A --> 1 : C --> 2 : G --> 3 : T --> 4)

    Args:
        text (str): DNA sequence

    Returns:
        array: Number array (representing the DNA sequence)
    """
    res_array = np.frombuffer(text.encode(), np.int8)
    return res_array


def strToBase(text: str) -> array:
    text = text.upper()
    res_array = np.empty(len(text))
    for i in range(len(text)):
        if text[i] == "A":
            res_array[i] = 1
        elif text[i] == "C":
            res_array[i] = 2
        elif text[i] == "G":
            res_array[i] = 3
        else:
            res_array[i] = 4
    return res_array


def asciiToStr(ascii_list: array) -> str:
    """Transform a Numpy array of int to a string
    TO DO : adapt to ACGT

    Args:
        ascii_list (array): Numpay array of int

    Returns:
        str: String representative of the initial array
    """
    res_txt = ascii_list.tobytes().decode("ascii")
    return res_txt


def makeP12(t_list: array) -> array:
    """Create P12 : Take all the n+1%3 and n+2%3 iif < N-2

    Args:
        t_list (array): array of number (corresponding to a DNA seq)

    Returns:
        array: array of all the n+1%3 and n+2%3 iif < N-2
    """
    lenght = len(t_list)
    p1 = np.arange(1, lenght - 2, 3)
    p2 = np.arange(2, lenght - 2, 3)
    return np.concatenate((p1, p2))


def makeTriplet(t_list: array, p12_list: array) -> array:
    """Create triplet with p12_array : Do of array of 3 from the T array respecting :
    [T[p12][i] ,T[p12][i] +1 ,T[p12][i] + 2] for each element of p12.

    Args:
        t_list (array): array of number (corresponding to a DNA seq)
        p12_list (array): array of all the n+1%3 and n+2%3 iif < N-2 for T

    Returns:
        array: Array of triplet create by p12 and T
    """
    r12 = np.zeros(len(p12_list), dtype=[("1", int), ("2", int), ("3", int)])
    for i, ind in enumerate(p12_list):
        r12[i] = (t_list[ind], t_list[ind + 1], t_list[ind + 2])
    return r12


def sortTriplet(r12: list[list[int]]) -> list[list[int]]:
    return quickSort(r12)


def makeDict(p12: array, r12: array, r12s: array, index_list: array, tp: array) -> dict:
    """Make a dict to keep all data needed for this algorithm
    TO DO : Deep only data needed

    Args:
        p12 (array): p12 array (Index of T kept with n+1 and n+2 % 3)
        r12 (array): Triplet array no sorted
        r12s (array): Triplet array sorted
        index_list (array): Array of index of each triplet before the sort
        tp (array): Array of index for each triplet after the sort, when it was before

    Returns:
        dict: Dict contening all data
    """
    index = np.take_along_axis(np.copy(p12), index_list, axis=0)
    i_dict = {
        "p12": p12,
        "r12": r12,
        "r12s": r12s,
        "index": index,
        "tp": tp,
    }
    return i_dict


def print_dict(iter_dict: dict) -> None:
    for key, vals in iter_dict.items():
        print(f"The {key} is : {vals}")


def print_iters_dict(iters_dict: dict) -> None:
    for key, vals in iters_dict.items():
        print(f"For the iteration number --- {key} --- :\n")
        print_dict(vals)


def reccursive_sort_s11(init_seq, iter_dict, iter=0):
    T = np.concatenate((init_seq, np.zeros(3)))
    p12 = makeP12(T)
    r12 = makeTriplet(T, p12)
    r12s = np.copy(r12)
    index_list, tp = sortTriplet(r12s)
    d = makeDict(p12, r12, r12s, index_list, tp)
    iter_dict[str(iter)] = d

    if tp.shape == np.unique(tp).shape:  # step 2 ?
        print("on est la")
    else:
        i = iter + 1
        reccursive_sort_s11(iter_dict[str(iter)]["tp"], iter_dict, iter=i)


def step2():
    pass


def dc3(seq: str):
    iter_dict = {}
    T = strToAscii(seq)
    reccursive_sort_s11(T, iter_dict)
    print(iter_dict)
    print_iters_dict(iter_dict)


if __name__ == "__main__":
    # s = "abcabcacab"
    # iter_dict = {}

    # T = strToAscii(s)
    # print(T)
    # reccursive_sort_s11(T, iter_dict)
    # print(iter_dict)
    # print_iters_dict(iter_dict)
    print(strToBase("ATTAGCAGCC"))


# Important pas stocker tout les info, en théorie y'a juste order a stocker. (P0 a voir ptete besoin a la fin).
# Recoder le "radix" short algo
