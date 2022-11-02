# Radix sort in Python

import cProfile
from random import randint
import numpy as np


def countingSort(array, index_list, place):
    size = len(array)
    output = np.zeros(size, dtype=int)
    count = np.zeros(10, dtype=int)
    output_ind_list = np.zeros(size, dtype=int)
    print(output)
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

        output[count[index % 10] - 1] = array[i]
        output_ind_list[i] = count[index % 10] - 1
        print(array[i], output)
        print(count[index % 10] - 1, output_ind_list)
        count[index % 10] -= 1

        i -= 1

    for i in range(0, size):
        array[i] = output[i]
        index_list[i] = output_ind_list[i]

    return index_list


def radixSort(array):
    index_list = np.arange(len(array), dtype=int)

    max_element = np.amax(array)

    place = 1
    while max_element // place > 0:
        countingSort(array, index_list, place)
        place *= 10

    return index_list


if __name__ == "__main__":
    data = np.array([121, 432, 1, 3, 1, 5, 45, 788], dtype=int)
    CONSTANT_DATA = np.copy(data)
    print(data)
    a = radixSort(data)
    print(data, a)

    datatest = [randint(10, 99) for i in range(500000)]
    # cProfile.run("radixSort(datatest)")

    # cProfile.run("sorted(datatest)")
