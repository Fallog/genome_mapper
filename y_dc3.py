from array import array
import numpy as np
import tools


def make_p12(t_list: array) -> array:
    """Create P12 : Take all the n+1%3 and n+2%3 iif < N-2

    Args:
        t_list (array): array of number (corresponding to a DNA seq)

    Returns:
        array: array of all the n+1%3 and n+2%3 iif < N-2
    """
    lenght = t_list.size
    p1 = np.arange(1, lenght - 2, 3, dtype=int)
    p2 = np.arange(2, lenght - 2, 3, dtype=int)
    return np.concatenate((p1, p2))


def make_triplet(t_list: array, p12_list: array) -> array:
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


def make_r0(p0, tp, index):
    """Make R0 : list of list of 2 int. In the following form :
    r0[i][0] = T[p0[i]] with T = T or T' (tp) if not the first iteration
    r0[i][0] = index of p0[i]+1 in the list of index of R12 after sort

    Args:
        p0 (array): all index of T / T' i such that i%3=0
        tp (array): T or T' if not the first iteration
        index (array): list of index of R12 after sort

    Returns:
        array: The r0 list (follow rules above)
    """
    y = p0.size
    r0 = np.empty((y), dtype=[("1", int), ("2", int)])
    for i in range(y):
        r0[i][0] = int(tp[p0[i]])
        r0[i][1] = (
            np.where(index == p0[i] + 1)[0][0] + 1
        )  # TODO: Maybe find a better way
    return r0


def compute_dc3_variable(T, sorting_algorithm):
    """
    Compute all DC3 variable computed in the first step and return it

    Args:
        T (array) : Array of all the values

    Returns:
        p12 (array): p12 array (Index of T kept with n+1 and n+2 % 3)
        r12 (array): Triplet array no sorted
        r12s (array): Triplet array sorted
        index (array): list of index of R12 after sort
    """

    p12 = make_p12(T)
    r12 = make_triplet(T, p12)
    r12s = np.copy(r12)
    index_list = np.argsort(r12s, kind=sorting_algorithm, order=("1", "2", "3"))
    index = np.take_along_axis(p12, index_list, axis=0)

    r12s.sort(
        kind="quicksort", order=("1", "2", "3")
    )  # TODO: change the sort to optimize ?

    return p12, r12, r12s, index


def make_tp(r12s, r12):
    tp = np.empty(r12s.size, dtype=int)
    array_without_duplicate = np.unique(r12s)
    for i, elem in enumerate(r12):
        tp[i] = np.where(elem == array_without_duplicate)[0][0] + 1
    return tp


def reccursive_sort_s11(init_seq, sorting_algorithm, iter=0, print_var=False):
    """Do the reccusrive sort to eliminate all value in double

    Args:
        init_seq (array): T string
        iter (int, optional): The actual number of iteration. Defaults to 0.
        print_var(Bool, optional): If true, print var of DC3 (not usefull)

    return:
        array: Suffix array sorted
    """
    T = np.concatenate((init_seq, np.zeros(3)))
    p12, r12, r12s, index = compute_dc3_variable(T, sorting_algorithm)
    if r12.shape != np.unique(r12, axis=0).shape:
        tp = make_tp(r12s, r12)
        i = iter + 1
        index012 = reccursive_sort_s11(tp, sorting_algorithm, i)
        index = np.take_along_axis(p12, index012, axis=0)
    if __name__ == "__main__" and print_var:
        tools.printDc3Var(p12, r12, r12s, index, T, iter)  # To check
    if iter:
        p0 = np.arange(0, init_seq.size, 3)
        r0 = make_r0(p0, init_seq, index)
        index012 = step2(T, index, r0, p0)
        # Index012 is corresponding of index des index d'avant
        return index012
    else:
        p0 = np.arange(0, init_seq.size, 3)  # init_seq = T sans les 0 psq iter = 0
        r0 = make_r0(p0, init_seq, index)
        res = step2(T, index, r0, p0)
        return res


def step2(T, index_r12, r0, p0):
    """Step of merging to create a new index array updated for the next iteration

        Args:
            T (array): Array of int corresponding to the str sequence for the i iteration
    with the sentinels 0
            index_r12 (array): list of index of R12 after sort
            r0 (array): r0 array (check make_r0)
            p0 (array): all index of T / T' i such that i%3=0

        Returns:
            array: New index array updated for the next iteration
    """
    index_r0 = np.argsort(
        r0, kind="mergesort", order=("1", "2")
    )  # Give index after sort
    index0 = np.take_along_axis(
        p0, index_r0, axis=0
    )  # Take an array and new index and give new array
    index012 = merging_r0_r12(T, index0, index_r12)
    return index012


def get_smallest_index(T, val_12, val_0, index_r12):  # TODO: Maybe change to optimize
    """Algorithm to get the smallest value adapted to DC3 necessity
    Its needed to merge r0 and r12 when value are equals

    Args:
        T (array): Array of int corresponding to the str sequence for the i iteration
    with the sentinels 0
        val_12 (int): current value in p12
        val_0 (int): current value in p0
        index_r12 (array): list of index of R12 after sort

    Returns:
        tuple: (value,"0" if belong to p0, else "12" (such that correspond to p12))
    """
    T_12 = T[val_12]
    T_0 = T[val_0]
    if T_12 != T_0:
        if min(T_0, T_12) == T_0:
            return (val_0, "0")
        else:
            return (val_12, "12")
    else:
        if val_12 % 3 == 1 or val_0 % 3 == 1:
            if (
                min(
                    np.where(index_r12 == val_12 + 1)[0][0],
                    np.where(index_r12 == val_0 + 1)[0][0],
                )
                == np.where(index_r12 == val_0 + 1)[0][0]
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
                    min(
                        np.where(index_r12 == val_12 + 2)[0][0],
                        np.where(index_r12 == val_0 + 2)[0][0],
                    )
                    == np.where(index_r12 == val_0 + 2)[0][0]
                ):  # We check smallest value bc already ordered in r12
                    return (val_0, "0")
                else:
                    return (val_12, "12")


def merging_r0_r12(T, index_r0, index_r12):
    """Merge r0 and r12 index to update the good index for the next iteration

    Args:
        T (array): Array of int corresponding to the str sequence for the i iteration
    with the sentinels 0
        index_r0 (_type_): is r0 (check make_r0)
        index_r12 (array): list of index of R12 after sort

    Returns:
        array: index array updated for the next iteration
    """
    size12, size0 = index_r12.size, index_r0.size
    index_table_merged = np.empty(size12 + size0, dtype=int)
    i_r12, i_r0 = 0, 0
    while i_r12 < size12 and i_r0 < size0:
        val_12 = index_r12[i_r12]
        val_0 = index_r0[i_r0]
        val, val_type = get_smallest_index(T, val_12, val_0, index_r12)
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
    if T.size % 3 == 0:
        return index_table_merged
    else:
        return index_table_merged[1:]  # Eleminate centinel


def dc3(seq: str, is_seq_to_ascii: bool = False, sorting_algorithm: str = "mergesort"):
    if is_seq_to_ascii:
        T = tools.strToAscii(seq)
    else:
        T = tools.strToBase(seq)
    suffix_array = reccursive_sort_s11(T, sorting_algorithm)
    return suffix_array


if __name__ == "__main__":
    test_result_dic = {
        "abcabcacab": [8, 0, 3, 6, 9, 1, 4, 7, 2, 5],
        "ab": [0, 1],
        "abaaba": [5, 2, 3, 0, 4, 1],
    }
    for i, (key, val) in enumerate(test_result_dic.items()):
        test = dc3(key, is_seq_to_ascii=True)
        print(
            f"""\n--- TEST {i+1} ---
Test with n%3=1 char string. Suffix array obtened = {test}
Is equal to the theorical result ? {(test == val).all()} !"""
        )
