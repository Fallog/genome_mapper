# import numpy   ## A convertir en numpy
from array import array
import numpy as np
import tools


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


def reccursive_sort_s11(init_seq, iter=0):
    T = np.concatenate((init_seq, np.zeros(3)))
    p12, r12, r12s, index, tp = computeDc3Variable(T)

    if tp.shape == np.unique(tp).shape:  # step 2 ?
        print("on est la")
    else:
        i = iter + 1
        reccursive_sort_s11(tp, i)
    print(tp, iter)  # STILL WORKING !


def computeDc3Variable(T):
    """
    Compute all DC3 variable computed in the first step and return it

    Args:
        T (array) : Array of all the values

    Returns:
        p12 (array): p12 array (Index of T kept with n+1 and n+2 % 3)
        r12 (array): Triplet array no sorted
        r12s (array): Triplet array sorted
        index (array): TODO: FIND WAY TO EXPLAIN THIS
        tp (array): Array of index for each triplet after the sort, when it was before (T\')
    """

    p12 = makeP12(T)
    r12 = makeTriplet(T, p12)
    r12s = np.copy(r12)
    index_list = np.argsort(r12s, kind="mergesort", order=("1", "2", "3"))
    index = np.take_along_axis(np.copy(p12), index_list, axis=0)

    r12s.sort(kind="quicksort", order=("1", "2", "3"))  # SORT HERE
    tp = np.empty(r12s.size)
    array_without_duplicate = np.unique(r12s)
    for i, elem in enumerate(r12):
        tp[i] = np.where(elem == array_without_duplicate)[0][0] + 1

    return p12, r12, r12s, index, tp


def step2():
    pass


def dc3(seq: str):  # TODO: Changer recursive_sort_s11 en DC3
    T = tools.strToAscii(seq)
    reccursive_sort_s11(T)
    # print(iter_dict)
    # print_iters_dict(iter_dict)


if __name__ == "__main__":
    # s = "abcabcacab"
    # iter_dict = {}

    # T = strToAscii(s)
    # print(T)
    # reccursive_sort_s11(T, iter_dict)
    # print(iter_dict)
    # print_iters_dict(iter_dict)
    iter_dict = {}
    iter_dict = dc3("abcabcacab")


# Important pas stocker tout les info, en théorie y'a juste order a stocker. (P0 a voir ptete besoin a la fin).
# Recoder le "radix" short algo
