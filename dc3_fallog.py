import cProfile
import random
from radix_sort_tuple_fallog import radixSort


def dnaSeqToIntTable(dna_seq: str) -> list[int]:
    """Convert a string of character into a list of integers representing each character.
    Because the strings are DNA sequence, the conversion is made with a dictionary mapping each nucleotide with a number.
    A -> 1 ; C -> 2 ; G -> 3 ; T -> 4
    The chosen number respect the order of this 4 letters in the alphabet.

    Args:
        dna_seq (str): sequence of A, C, G, T

    Returns:
        list[int]: representation of dna_seq with an array of integers
    """
    dna_dict = {"A": 1, "C": 2, "G": 3, "T": 4}
    int_array = [0] * len(dna_seq)
    for index, nucleotide in enumerate(dna_seq):
        int_array[index] = dna_dict[nucleotide]
    return int_array


def appendSentinelNumbers(input_list: list[int]) -> list:
    """Returns input_list with three 0 appended at the end of the list.

    Args:
        input_list (list[int]): a list of integer

    Returns:
        list: input_list but concatenated with a list of three 0
    """
    return input_list + [0, 0, 0]


def getPositionTable(input_list: list[int]) -> list[int]:
    """Returns the concatenation of an array containing indexes i of input_list that satisfy i % 3 = 1 and
    another array containing indexes satisfying i % 3 = 2. Also i = 2 < len(input_list)

    Args:
        input_list (list[int]): an array of integers

    Returns:
        list[int]: concatenation of 2 arrays, one containing indexes i of input_list such that i % 3 = 1 and another
        such that i % 3 = 2.
    """
    size = len(input_list)
    P1 = [i for i in range(size) if i % 3 == 1 and i + 2 < size]
    P2 = [i for i in range(size) if i % 3 == 2 and i + 2 < size]
    return P1 + P2


def getTripletList(
    input_list: list[int], position_table: list[int]
) -> list[tuple[tuple, int]]:
    pos_size = len(position_table)
    triplet_list = [0] * pos_size
    for i in range(pos_size):
        position = position_table[i]
        triplet_list[i] = (
            (
                input_list[position],
                input_list[position + 1],
                input_list[position + 2],
            ),
            1,
        )
    return triplet_list


if __name__ == "__main__":
    test_seq = "ATTAGCAGCC"
    print("Length of the sequence: ", len(test_seq))
    seq_int = dnaSeqToIntTable(test_seq)
    print(seq_int == [1, 4, 4, 1, 3, 2, 1, 3, 2, 2])  # OK

    seq_int = appendSentinelNumbers(seq_int)
    print(seq_int == [1, 4, 4, 1, 3, 2, 1, 3, 2, 2, 0, 0, 0])  # OK

    pos_table = getPositionTable(seq_int)
    print(pos_table == [1, 4, 7, 10, 2, 5, 8])  # OK

    triplet_list = getTripletList(seq_int, pos_table)
    print(
        triplet_list
        == [
            ((4, 4, 1), 1),
            ((3, 2, 1), 1),
            ((3, 2, 2), 1),
            ((0, 0, 0), 1),
            ((4, 1, 3), 1),
            ((2, 1, 3), 1),
            ((2, 2, 0), 1),
        ]
    )  # OK
    # sorted_triplet_list = radixSort(triplet_list)
