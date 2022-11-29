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
    lenght = t_list.size
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
    r12 = np.zeros(p12_list.size, dtype=[("1", int), ("2", int), ("3", int)])
    for i, ind in enumerate(p12_list):
        r12[i] = (t_list[ind], t_list[ind + 1], t_list[ind + 2])
    return r12


def makeR0(p0, tp):
    y = p0.size
    r0 = np.empty((y), dtype=[("1", int), ("2", int)])
    for i in range(y):  # TODO: TO UPGRADE
        r0[i][0] = int(tp[p0[i]])
        r0[i][1] = p0[i]  # NOT SURE
    return r0


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
    index = np.take_along_axis(p12, index_list, axis=0)

    r12s.sort(kind="quicksort", order=("1", "2", "3"))  # SORT HERE
    tp = np.empty(r12s.size)
    array_without_duplicate = np.unique(r12s)
    for i, elem in enumerate(r12):
        tp[i] = np.where(elem == array_without_duplicate)[0][0] + 1

    return p12, r12, r12s, index, tp


def reccursive_sort_s11(init_seq, last_tp=None, iter=0):
    """Do the reccusrive sort to eliminate all value in double

    Args:
        init_seq (array): T string
        last_tp (_type_): Last T\'.
        iter (int, optional): The actual number of iteration. Defaults to 0.
    """
    T = np.concatenate((init_seq, np.zeros(3)))
    if last_tp is None:
        last_tp = T
    p12, r12, r12s, index, tp = computeDc3Variable(T)
    if tp.shape == np.unique(tp).shape:  # step 2 ?
        print("on est la")
    else:
        i = iter + 1
        reccursive_sort_s11(tp, tp, i)
    tools.printDc3Var(p12, r12, r12s, index, tp, iter)  ## Pour vérifier

    if iter:
        step2(T, index, init_seq, iter)


def step2(T, index_r12, tp, iter):
    p0 = np.arange(0, tp.size, 3)
    r0 = makeR0(p0, tp)
    index_r0 = np.argsort(
        r0, kind="mergesort", order=("1", "2")
    )  # Give index after sort
    index0 = np.take_along_axis(
        p0, index_r0, axis=0
    )  # Take an array and new index and give new array
    index0123 = merging_r0_r12(T, index0, index_r12)
    print(index0123)


def get_smallest_index(T, val_12, val_0):  # Maybe too long
    T_12 = T[val_12]
    T_0 = T[val_0]
    if T_12 != T_0:
        if min(T_0, T_12) == T_0:
            return (val_0, "0")
        else:
            return (val_12, "12")
    else:
        if val_12 % 3 == 1:
            if (
                min(val_12, val_0) == val_0
            ):  # We check smallest value bc already ordered in r12
                return (val_0, "0")
            else:
                return (val_12, "12")
        else:
            Tp_12 = T[val_12 + 1]
            Tp_0 = T[val_0 + 1]
            if Tp_12 != Tp_0:
                if min(Tp_0, Tp_12) == Tp_0:
                    return (val_0, "0")
                else:
                    return (val_12, "12")
            else:
                if (
                    min(val_12, val_0) == val_0
                ):  # We check smallest value bc already ordered in r12
                    return (val_0, "0")
                else:
                    return (val_12, "12")


def merging_r0_r12(T, index_r0, index_r12):
    size12, size0 = index_r12.size, index_r0.size
    index_table_merged = np.empty(size12 + size0)
    i_r12, i_r0 = 0, 0
    while i_r12 < size12 and i_r0 < size0:
        val_12 = index_r12[i_r12]
        val_0 = index_r0[i_r0]
        val, val_type = get_smallest_index(T, val_12, val_0)
        index_table_merged[i_r0 + i_r12] = val
        if val_type == "0":
            i_r0 += 1
        else:
            i_r12 += 1
    if i_r12 == size12:  # Case of r12 is finished
        while i_r0 < size0:
            index_table_merged[i_r12 + i_r0] = index_r0[i_r0]
            i_r0 += 1

    else:
        while i_r12 < size12:
            index_table_merged[i_r0 + i_r12] = index_r12[i_r12]
            i_r12 += 1
    return index_table_merged[1:]  # Eleminate centinel


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
