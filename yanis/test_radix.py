# Radix sort in Python

import cProfile
from random import randint
import numpy as np


def countingSort(array, place):
    size = len(array)
    output = [0] * size
    count = [0] * 10

    # Calculate count of elements
    for i in range(0, size):
        index = array[i] // place
        count[index % 10] += 1
    # print(count)
    # Calculate cumulative count
    for i in range(1, 10):
        count[i] += count[i - 1]
    # print(count)

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


# Main function to implement radix sort
def radixSort(array):
    data2 = array[:]

    max_element = max(array)

    place = 1
    while max_element // place > 0:
        countingSort(array, place)
        place *= 10

    list_without_duplicate = list(set(array))
    init_index = [list_without_duplicate.index(i) for i in data2]

    return init_index


data = [121, 432, 1, 3, 1, 5, 45, 788]
print(data)
a = radixSort(data)
print(data, a)


datatest = [randint(90, 99) for i in range(500000)]
cProfile.run("radixSort(datatest)")
