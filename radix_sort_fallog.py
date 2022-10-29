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


def generateCountArray(digits_list: list[int], cumulative: bool = False) -> list[int]:
    """Computes an array containing the number of occurences of every digits in digits_list.


    Args:
        digits_list (list[int]): list of integers of only one digit long (from 0 to 9)
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


def createDigitList(int_list: list[int], target_digit: int) -> list[int]:
    """Extract the target_digit in a list of integers

    Args:
        int_list (list[int]): every integers which digits are extracted
        target_digit (int): position of the extracted digit, 1 refering to the unit digit in the integer,
        2 for the tens, 3 for hundreds etc.

    Returns:
        list[int]: every extracted digits
    """
    digits_list = []
    for integer in int_list:
        if not extractDigit(integer, target_digit) < 1:
            digits_list.append(extractDigit(integer, target_digit))
        else:
            digits_list.append(0)
    return digits_list


def countingSort(int_list: list[int], target_digit: int) -> list[int]:
    """Sort an array of integers in ascending order based only on the target_digit
    Args:
        int_list (list[int]): array of integers

    Returns:
        list[int]: sorted version of int_list, in ascending order
    """
    count_array = generateCountArray(
        createDigitList(int_list, target_digit), cumulative=True
    )
    sorted_array = [0] * len(int_list)
    for integer in int_list[::-1]:  # reverse int_list
        new_index = count_array[extractDigit(integer, target_digit)]
        sorted_array[new_index - 1] = integer
        count_array[extractDigit(integer, target_digit)] -= 1
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
    print(getDigitsNumber(9) == 1)  # OK
    print(getDigitsNumber(53) == 2)  # OK
    print(getDigitsNumber(54136374251) == 11)  # OK

    print(
        generateCountArray([1, 1, 1, 6, 5, 7, 2, 6]) == [0, 3, 1, 0, 0, 1, 2, 1]
    )  # OK
    print(
        generateCountArray([1, 1, 1, 6, 5, 7, 2, 6], cumulative=True)
        == [0, 3, 4, 4, 4, 5, 7, 8]
    )  # OK

    print(extractDigit(7, 2) == 0)  # OK
    print(extractDigit(1452612, 1) == 2)  # OK

    test_list = [121, 432, 564, 23, 1, 45, 788]
    print(createDigitList(test_list, 1) == [1, 2, 4, 3, 1, 5, 8])  # OK
    print(createDigitList(test_list, 2) == [2, 3, 6, 2, 0, 4, 8])  # OK
    print(createDigitList(test_list, 3) == [1, 4, 5, 0, 0, 0, 7])  # OK

    digit_test_list = [1, 2, 4, 3, 1, 5, 8]
    countingSort(test_list, 1)

    print(countingSort(test_list, 1) == [121, 1, 432, 23, 564, 45, 788])  # OK
    print(countingSort(test_list, 2) == [1, 121, 23, 432, 45, 564, 788])  # OK
    print(countingSort(test_list, 3) == [23, 1, 45, 121, 432, 564, 788])  # OK

    print(radixSort(test_list) == [1, 23, 45, 121, 432, 564, 788])

    large_test = [random.randint(0, 5000) for i in range(50000)]
    cProfile.run("radixSort(large_test)")
