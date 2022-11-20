from array import array
import cProfile
from random import randint, seed
import numpy as np


def quickSort(x):
    """
    Return index list and T\'. The x array in parameter is sorted.
    """
    x_copy = np.copy(x)

    index_list = np.argsort(x, kind="mergesort", order=("1", "2", "3"))
    x.sort(kind="quicksort", order=("1", "2", "3"))
    tp = np.empty(x.size)
    array_without_duplicate = np.unique(x)
    for i, elem in enumerate(x_copy):
        tp[i] = np.where(elem == array_without_duplicate)[0][0] + 1

    return index_list, tp


if __name__ == "__main__":

    liste = np.array(
        [
            (5, 28, 1),
            (9, 111, 9),
            (9, 1, 5),
            (99, 1, 9),
            (2, 1, 89),
            (29, 19, 1),
        ],
        dtype=[("1", int), ("2", int), ("3", int)],
    )

    seed(12)
    # datatest = [(randint(1, 4), randint(1, 4), randint(1, 4)) for i in range(500000)]

    # datatest = np.array(datatest, dtype=[("1", int), ("2", int), ("3", int)])

    # print(datatest)

    print(liste)
    res = quickSort(liste)
    print(liste)
    print(res)
    # cProfile.run("quickSort(datatest)")
