###########
# IMPORTS #
###########
import math
import cProfile
import random
from IDTuple import IDTuple
from IDTupleList import IDTupleList


def getDigitsNumber(integer: int) -> int:
    """Computes the total number of digits of a given integer

    Args:
        integer (int): any number of type 'int'

    Returns:
        int: The total number of digits of
    """
    return math.ceil(math.log10(integer))


def findMaxinTupleList(tuple_list: tuple[int]) -> int:
    greatest_int = 0
    for tuple in tuple_list:
        max_element = max(tuple)
        if max_element > greatest_int:
            greatest_int = max_element
    return greatest_int


def generateCountArray(
    id_tuple_l: list[IDTuple], cumulative: bool = False
) -> list[int]:
    """Computes an array containing the number of occurences of every digits in digits_list.


    Args:
        digits_list (list[tuple]): list of integers of only one digit long (from 0 to 9)
        cumulative (bool): True if the count is done cumulatively ([1,0,2,1,4] -> [1,3,4,4,5]), False either

    Returns:
        list[int]: list of the count of the digits present in digits_lits
    """
    count_array = list(range(findMaxinTupleList(id_tuple_l) + 1))
    if cumulative:
        total_digits = 0
        for digit in count_array:
            total_digits += id_tuple_l.count(digit)
            count_array[digit] = total_digits
    else:
        for digit in count_array:
            count_array[digit] = id_tuple_l.count(digit)
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


def createLinkedDigitList(tuple_list: list[int], target_digit: int) -> list[int]:
    """Extract the target_digit in a list of integers

    Args:
        int_list (list[int]): every integers which digits are extracted
        target_digit (int): position of the extracted digit, 1 refering to the unit digit in the integer,
        2 for the tens, 3 for hundreds etc.

    Returns:
        list[int]: every extracted digits
    """
    output_list = [0] * len(tuple_list)
    for index, tuple in enumerate(tuple_list):
        digit = extractDigit(tuple[0], target_digit)
        output_list[index] = (digit, tuple[1])
    return output_list


def extractIntsFromTuples(tuple_list: list[tuple[int]], wanted_int: int) -> list[int]:
    size = len(tuple_list)
    wanted_ints_list = [0] * size
    for index, tuple in enumerate(tuple_list):
        wanted_ints_list[index] = tuple[(wanted_int - 1) % len(tuple)]

    return wanted_ints_list


def countingSort(int_list: list[int], target_digit: int):
    """Sort an array of integers in ascending order based only on the target_digit
    Args:
        int_list (list[int]): array of integers

    Returns:
        list[int]: sorted version of int_list, in ascending order
    """
    count_array = generateCountArray(
        createLinkedDigitList(int_list, target_digit), cumulative=True
    )
    size = len(int_list)
    sorted_array = [0] * size
    i = size - 1
    while i >= 0:
        integer = int_list[i]
        digit = extractDigit(integer, target_digit)
        new_index = count_array[digit]
        sorted_array[new_index - 1] = integer
        count_array[digit] -= 1
        i -= 1
    return sorted_array


def radixSort(int_list: list[int]) -> None:
    """Sort an array of integers in ascending order.

    Args:
        int_list (list[int]): array of integers

    Returns:
        _type_: sorted version of int_list, in ascending order
    """
    nb_loops = getDigitsNumber(max(int_list))
    for target_digit in range(1, nb_loops + 1):
        sorted_list = countingSort(int_list, target_digit)
        int_list = sorted_list
    return sorted_list


if __name__ == "__main__":
    ##############
    # UNIT TESTS #
    ##############
    large_test = [random.randint(0, 5000) for i in range(500000)]
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
    print(extractIntsFromTuples(tuple_list, 1) == [98, 98, 99, 0, 99, 99, 97])
    f_id_tuple = IDTupleList(extractIntsFromTuples(tuple_list, 1))
    s_id_tuple = IDTupleList(extractIntsFromTuples(tuple_list, 2))
    t_id_tuple = IDTupleList(extractIntsFromTuples(tuple_list, 3))
