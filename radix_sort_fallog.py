###########
# IMPORTS #
###########
import math


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
    count_array = list(range(10))
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
        target_digit (int): an integer > 1 -> position of the extracted digit, 1 refering to the rightest digit in the integer

    Returns:
        int: the digit to be extracted
    """
    return int(str(integer // 10 ** (target_digit - 1))[-1])


def createDigitList(int_list: list[int], target_digit: int) -> list[int]:
    """Extract the target_digit in a list of integers

    Args:
        int_list (list[int]): every integers which digits are extracted
        target_digit (int): position of the extracted digit, 1 refering to the rightest digit in the integer

    Returns:
        list[int]: every extracted digits
    """
    digits_list = []
    for integer in int_list:
        if not integer // 10 ** (target_digit - 1) < 1:
            digits_list.append(extractDigit(integer, target_digit))
    return digits_list


def radixSort(int_list: list[int]) -> None:
    print(f"Starting list {int_list}")
    sorted_list = list(range(len(int_list)))
    nb_loops = getDigitsNumber(max(int_list))
    for i in range(1, nb_loops + 1):
        digits_list = createDigitList(int_list, i)
        count_array = generateCountArray(digits_list, cumulative=True)
        print(f"Digits list {digits_list}")
        print(f"Count array {count_array}")
        for j in range(1, len(digits_list) + 1):
            print(int_list)
            actual_digit = digits_list[-j]
            if not actual_digit == 0:
                nb_same_digit = digits_list[digits_list.index(actual_digit) :].count(
                    actual_digit
                )
                value_in_count_array = count_array[actual_digit]
                print(
                    f"Digit {actual_digit} | count_array[{actual_digit}] {value_in_count_array}"
                )
                count_array[actual_digit] = value_in_count_array - nb_same_digit
                sorted_list[value_in_count_array - nb_same_digit] = int_list[-j]
                print(
                    f"Index {value_in_count_array - nb_same_digit} is set to {int_list[-j]}"
                )
        print(f"List after loop {sorted_list}")
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
        generateCountArray([1, 1, 1, 6, 5, 7, 2, 6]) == [0, 3, 1, 0, 0, 1, 2, 1, 0, 0]
    )  # OK
    print(
        generateCountArray([1, 1, 1, 6, 5, 7, 2, 6], cumulative=True)
        == [0, 3, 4, 4, 4, 5, 7, 8, 8, 8]
    )  # OK

    print(extractDigit(7, 2) == 0)  # OK
    print(extractDigit(1452612, 1) == 2)  # OK

    test_list = [56, 32, 9, 21, 991, 484, 17]
    print(createDigitList(test_list, 1) == [6, 2, 9, 1, 1, 4, 7])  # OK
    print(createDigitList(test_list, 2) == [5, 3, 2, 9, 8, 1])  # OK
    print(createDigitList(test_list, 3) == [9, 4])  # OK

    print(radixSort(test_list))
