import cProfile
import random
from radix_sort_fallog import radixSort


def strToAsciiTable(input_string: str) -> list[int]:
    """Convert a string of character into a list of integers representing each character in its ascii code

    Args:
        input_string (str): sequence of characters

    Returns:
        list[int]: list of the character's ascii code
    """
    ascii_table = [0] * len(input_string)
    for index, char in enumerate(input_string):
        ascii_table[index] = ord(char)
    return ascii_table


def appendSentinelNumbers(input_list: list[int]) -> list:
    return input_list + [0, 0, 0]


def getPositionTable(input_list: list[int]) -> list[int]:
    size = len(input_list)
    P1 = [i for i in range(size) if i % 3 == 1 and i + 2 < size]
    P2 = [i for i in range(size) if i % 3 == 2 and i + 2 < size]
    return P1 + P2


def getTripletList(input_list: list[int], position_table: list[int]) -> list[tuple]:
    pos_size = len(position_table)
    triplet_list = [0] * pos_size
    for i in range(pos_size):
        position = position_table[i]
        triplet_list[i] = (
            input_list[position],
            input_list[position + 1],
            input_list[position + 2],
        )
    return triplet_list


if __name__ == "__main__":
    # A -> 65, C -> 67, G -> 71, T -> 84  DNA bases in ASCII code
    test_seq = "abcabcacab"
    seq_ascii = strToAsciiTable(test_seq)
    print(seq_ascii == [97, 98, 99, 97, 98, 99, 97, 99, 97, 98])  # OK

    seq_ascii = appendSentinelNumbers(seq_ascii)
    print(seq_ascii == [97, 98, 99, 97, 98, 99, 97, 99, 97, 98, 0, 0, 0])  # OK

    pos_table = getPositionTable(seq_ascii)
    print(pos_table == [1, 4, 7, 10, 2, 5, 8])  # OK

    triplet_list = getTripletList(seq_ascii, pos_table)
    print(
        triplet_list
        == [
            (98, 99, 97),
            (98, 99, 97),
            (99, 97, 98),
            (0, 0, 0),
            (99, 97, 98),
            (99, 97, 99),
            (97, 98, 0),
        ]
    )  # OK
