# Radix sort in Python

import cProfile
from random import randint
import numpy as np


def countingSort(array, place):
    size = len(array)
    output = np.empty(size, dtype=int)
    count = np.zeros(10, dtype=int)

    # Calculate count of elements
    for i in range(0, size):
        index = array[i] // place
        count[index % 10] += 1

    # Calculate cumulative count
    for i in range(1, 10):
        count[i] += count[i - 1]

    # Place the elements in sorted order
    i = size - 1
    while i >= 0:
        index = array[i] // place
        # print(index)
        output[count[index % 10] - 1] = array[i]
        count[index % 10] -= 1
        i -= 1

    for i in range(0, size):
        array[i] = output[i]


def radixSort(array):
    data2 = np.copy(array)

    max_element = np.amax(array)

    place = 1
    while max_element // place > 0:
        countingSort(array, place)
        place *= 10

    list_without_duplicate = np.unique(array)

    init_index = np.empty_like(array)
    for i, elem in enumerate(data2):
        init_index[i] = np.where(list_without_duplicate == elem)[0][0]

    return init_index


if __name__ == "__main__":
    data = np.array([121, 432, 1, 3, 1, 5, 45, 788], dtype=int)
    print(data)
    a = radixSort(data)
    print(data, a)

    datatest = [randint(1000000, 11000000) for i in range(50000)]
    cProfile.run("radixSort(datatest)")

    cProfile.run("sorted(datatest)")
