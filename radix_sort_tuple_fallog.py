###########
# IMPORTS #
###########
import math
import cProfile
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
    tuple_list: list[tuple[tuple, int]], tuple_index: int, target_digit: int
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
        digit = extractDigit(tuple[0][tuple_index], target_digit)
        if not digit < 1:
            digits_list[index] = digit
    return digits_list


def findMaxElementTupleList(tuple_list: list[tuple[tuple, int]]):
    max_elem = 0
    for tuple in tuple_list:
        max_in_tuple = max(tuple[0])
        if max_in_tuple > max_elem:
            max_elem = max_in_tuple
    return max_elem


def countingSort(
    tuple_list: list[tuple[tuple, int]], target_digit: int
) -> list[tuple[tuple, int]]:
    """Sort an array of integers in ascending order based only on the target_digit
    Args:
        int_list (list[int]): array of integers

    Returns:
        list[int]: sorted version of int_list, in ascending order
    """
    tuple_size = len(tuple_list[0][0])
    for tuple_index in range(tuple_size):
        count_array = generateCountArray(
            createDigitList(tuple_list, tuple_index, target_digit), cumulative=True
        )
        size = len(tuple_list)
        sorted_array = [((0, 0, 0), 0)] * size
        i = size - 1
        while i >= 0:
            actual_tuple = tuple_list[i][0]
            integer = actual_tuple[tuple_index]
            digit = extractDigit(integer, target_digit)
            new_index = count_array[digit]
            sorted_array[new_index - 1] = (actual_tuple, new_index)
            count_array[digit] -= 1
            i -= 1
        return sorted_array


def radixSort(tuple_list: list[tuple[tuple, int]]) -> list[tuple[tuple, int]]:
    """Sort an array of integers in ascending order.

    Args:
        int_list (list[int]): array of integers

    Returns:
        _type_: sorted version of int_list, in ascending order
    """
    nb_loops = getDigitsNumber(findMaxElementTupleList(tuple_list))
    for target_digit in range(1, nb_loops + 1):
        sorted_list = countingSort(tuple_list, target_digit)
        tuple_list = sorted_list
    orderTupleList(tuple_list)
    return sorted_list


def orderTupleList(tuple_list: list[tuple[tuple, int]]):
    for i in range(1, len(tuple_list)):
        parsed_tuple = tuple_list[i][0]
        previous_tuple = tuple_list[i - 1][0]
        if parsed_tuple == previous_tuple:
            tuple_list[i] = (parsed_tuple, tuple_list[i - 1][1])
        else:
            tuple_list[i] = (parsed_tuple, tuple_list[i - 1][1] + 1)


if __name__ == "__main__":
    ##############
    # UNIT TESTS #
    ##############
    tuple_list = [
        ((98, 99, 97), 1),
        ((98, 99, 97), 1),
        ((99, 97, 98), 1),
        ((0, 0, 0), 1),
        ((99, 97, 98), 1),
        ((99, 97, 99), 1),
        ((97, 98, 0), 1),
    ]
    print(findMaxElementTupleList(tuple_list) == 99)  # OK
    print(createDigitList(tuple_list, 0, 1) == [8, 8, 9, 0, 9, 9, 7])  # OK
    print(
        generateCountArray([8, 8, 9, 0, 9, 9, 7], cumulative=True)
        == [1, 1, 1, 1, 1, 1, 1, 2, 4, 7]
    )  # OK
    print(
        countingSort(countingSort(tuple_list, 1), 2)
        == [
            ((0, 0, 0), 1),
            ((97, 98, 0), 2),
            ((98, 99, 97), 3),
            ((98, 99, 97), 4),
            ((99, 97, 98), 5),
            ((99, 97, 98), 6),
            ((99, 97, 99), 7),
        ]
    )  # OK
    print(
        radixSort(tuple_list)
        == [
            ((0, 0, 0), 1),
            ((97, 98, 0), 2),
            ((98, 99, 97), 3),
            ((98, 99, 97), 3),
            ((99, 97, 98), 4),
            ((99, 97, 98), 4),
            ((99, 97, 99), 5),
        ]
    )  # OK

    random.seed(14)
    large_tuple_list = [
        ((random.randint(0, 5000), random.randint(0, 5000), random.randint(0, 5000)), 0)
        for i in range(5000000)
    ]
    cProfile.run("radixSort(large_tuple_list)")
