from array import array
import cProfile
from random import randint
import numpy as np


def newSort(x):
    x_copy = np.copy(x)
    x.sort()
    list_index = np.empty(x.size)
    array_without_duplicate = np.unique(x)
    for i, elem in enumerate(x_copy):
        list_index[i] = np.where(elem == array_without_duplicate)[0][0]
    return list_index


liste = np.array(
    [(541, 25, 1), (1, 57, 126), (80, 80, 80), (70, 80, 90), (1, 57, 126), (2, 1, 1)],
    dtype=[("1", int), ("2", int), ("3", int)],
)


datatest = [randint(1000000, 11000000) for i in range(50000)]
datatest = np.array(datatest)
newSort(liste)

print(liste)

cProfile.run("newSort(datatest)")
