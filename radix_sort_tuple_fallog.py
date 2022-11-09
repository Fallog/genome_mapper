###########
# IMPORTS #
###########
import math

# import cProfile
import random


def getDigitsNumber(integer: int) -> int:
    """Computes the total number of digits of a given integer

    Args:
        integer (int): any number of type 'int'

    Returns:
        int: The total number of digits of
    """
    return math.ceil(math.log10(integer))


def generateCountArray(digits_list: list[tuple], cumulative: bool = False) -> list[int]:
    """Computes an array containing the number of occurences of every digits in digits_list.


    Args:
        digits_list (list[tuple]): list of integers of only one digit long (from 0 to 9)
        cumulative (bool): True if the count is done cumulatively ([1,0,2,1,4] -> [1,3,4,4,5]), False either

    Returns:
        list[int]: list of the count of the digits present in digits_lits
    """
    count_array = list(range(max(digits_list) + 1))
    if cumulative:
        total_digits = 0
        for digit in count_array:
            total_digits += digits_list.count(digit)
            count_array[digit] = total_digits
    else:
        for digit in count_array:
            count_array[digit] = digits_list.count(digit)
    return count_array


def extractDigit(integer: int, target_digit: int) -> int:
    """Extracts the target_digit in integer

    Args:
        integer (int): number of type 'int' which digit is extracted
        target_digit (int): an integer > 1 -> position of the extracted digit, 1 refering to the unit digit in the integer,
        2 for the tens, 3 for hundreds etc.

    Returns:
        int: the digit to be extracted
    """
    return (integer // 10 ** (target_digit - 1)) % 10


def createDigitList(
    tuple_list: list[tuple[int]], tuple_index: int, target_digit: int
) -> list[int]:
    """Extract the target_digit in a list of integers

    Args:
        tuple_list (list[tuple[int]]): list of triplets
        target_digit (int): position of the extracted digit, 1 refering to the unit digit in the integer,
        2 for the tens, 3 for hundreds etc.

    Returns:
        list[int]: every extracted digits
    """
    digits_list = [0] * len(tuple_list)
    for index, tuple in enumerate(tuple_list):
        digit = extractDigit(tuple[tuple_index], target_digit)
        if not digit < 1:
            digits_list[index] = digit
    return digits_list


def findMaxElementTupleList(tuple_list: list[tuple[int]]):
    max_elem = 0
    for tuple in tuple_list:
        max_elem = max(max_elem, max(tuple))
    return max_elem


def countingSort(tuple_list: list[tuple[int]], target_digit: int):
    """Sort an array of integers in ascending order based only on the target_digit
    Args:
        int_list (list[int]): array of integers

    Returns:
        list[int]: sorted version of int_list, in ascending order
    """
    tuple_size = len(tuple_list[0])
    for tuple_index in range(tuple_size):
        count_array = generateCountArray(
            createDigitList(tuple_list, tuple_index, target_digit), cumulative=True
        )
        size = len(tuple_list)
        sorted_array = [0] * size
        i = size - 1
        while i >= 0:
            integer = tuple_list[i][tuple_index]
            digit = extractDigit(integer, target_digit)
            new_index = count_array[digit]
            sorted_array[new_index - 1] = tuple_list[i]
            count_array[digit] -= 1
            i -= 1
        return sorted_array


def radixSort(tuple_list: list[tuple[int]]) -> None:
    """Sort an array of integers in ascending order.

    Args:
        int_list (list[int]): array of integers

    Returns:
        _type_: sorted version of int_list, in ascending order
    """
    nb_loops = getDigitsNumber(max(tuple_list))
    for target_digit in range(1, nb_loops + 1):
        sorted_list = countingSort(tuple_list, target_digit)
        tuple_list = sorted_list
    return sorted_list


if __name__ == "__main__":
    ##############
    # UNIT TESTS #
    ##############
    large_tuple_list = [random.randint(0, 5000) for i in range(500000)]
    # cProfile.run("radixSort(large_test)")
    tuple_list = [
        (98, 99, 97),
        (98, 99, 97),
        (99, 97, 98),
        (0, 0, 0),
        (99, 97, 98),
        (99, 97, 99),
        (97, 98, 0),
    ]
    print(findMaxElementTupleList([(0, 66, 89), (51, 23, 18)]) == 89)  # OK
    print(createDigitList(tuple_list, 0, 1) == [8, 8, 9, 0, 9, 9, 7])  # OK
    print(
        generateCountArray([8, 8, 9, 0, 9, 9, 7], cumulative=True)
        == [1, 1, 1, 1, 1, 1, 1, 2, 4, 7]
    )  # OK
    print(
        countingSort(countingSort(tuple_list, 1), 2)
        == [
            (0, 0, 0),
            (97, 98, 0),
            (98, 99, 97),
            (98, 99, 97),
            (99, 97, 98),
            (99, 97, 98),
            (99, 97, 99),
        ]
    )  # OK
