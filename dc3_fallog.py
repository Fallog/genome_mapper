# import cProfile
# import random
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
    another array containing indexes satisfying i % 3 = 2. Also i + 2 < len(input_list).

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


def getTripletListIndexed(
    input_list: list[int], position_table: list[int]
) -> list[tuple[tuple[int, int, int], int]]:
    """Returns a list of tuple made of a triplet plus its index in the position table.
    The triplets are built using the position table and an array of integers.
    For each element p in the position table, we form a tuple looking like this :
    ((input_list[p], input_list[p + 1], input_list[p + 2]), p)

    Args:
        input_list (list[int]): an array of integers
        position_table (list[int]): array obtained by the getPositionTable function

    Returns:
        list[tuple[tuple[int, int, int], int]]: list of tuple where a tuple consists in a triplet and
    an integer representing the index from which the triplet is built
    """
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
            position,
        )
    return triplet_list


def getTripletList(
    index_triplet_list: list[tuple[tuple, int]]
) -> list[tuple[int, int, int]]:
    """Returns a list of triplets, extracted from a list of tuple where a tuple consists in a triplet and
    an integer representing the order of the triplet

    Args:
        index_triplet_list (list[tuple[tuple, int]]): list of tuple where a tuple consists in a triplet and
    an integer representing the index from which the triplet was built in the first place

    Returns:
        list[tuple[int, int, int]]: list made of the triplets of ordered_triplet_list without the index
    """
    return [index_triplet_list[i][0] for i in range(len(index_triplet_list))]


def getIndexTable(triplet_list: list[tuple[tuple[int, int, int], int]]) -> list[int]:
    """Returns an array of integers built with the int inside the tuple of triplet_list.
    triplet_list is an array of ((int, int, int), int) and here it extracts the last int and put it in the returned list.

    Args:
        triplet_list (list[tuple[tuple[int, int, int], int]]): a tuple consists in a triplet and an integer representing
        the index from which the triplet was built in the first place

    Returns:
        list[int]: indexes of the triplet in an sorted or unsorted triplet_list
    """
    return [triplet_list[i][1] for i in range(len(triplet_list))]


def getOrderTable(
    sorted_tuple_list: list[tuple[tuple[int, int, int], int]]
) -> list[int]:
    """Returns the sorted position of each triplet in sorted_tuple_list.
    If two triplet are identical, their position is set equal.

    Args:
        sorted_tuple_list (list[tuple[tuple[int, int, int], int]]): a tuple consists in a triplet and an integer
        representing the index from which the triplet was built in the first place

    Returns:
        list[int]: sorted positions of the triplets
    """
    size = len(tuple_list)
    order_table = [1] * size

    for i in range(1, size):
        parsed_tuple = sorted_tuple_list[i][0]
        previous_tuple = sorted_tuple_list[i - 1][0]
        if parsed_tuple == previous_tuple:
            order_table[i] = order_table[i - 1]
        else:
            order_table[i] = order_table[i - 1] + 1
    return order_table


def createNextIntList(
    index_table: list[int], index_table_sort: list[int], order_table: list[int]
) -> list[int]:
    """Returns an array made of the index + 1 in index_table_sort of each element of index_table.

    Args:
        index_table (list[int]): indexes of the triplet in an unsorted triplet_list
        index_table_sort (list[int]): indexes of the triplet in a sorted triplet_list
        order_table (list[int]): orders of the triplet in a sorted triplet_list

    Returns:
        list[int]:
    """
    size = len(index_table)
    next_int_list = [0] * size
    for i in range(size):
        next_int_list[i] = order_table[index_table_sort.index(index_table[i])]
    return next_int_list


def DC3(int_sequence: list[int]):
    pass


if __name__ == "__main__":
    test_seq = "ATTAGCAGCC"
    print("Length of the sequence: ", len(test_seq))
    seq_int = dnaSeqToIntTable(test_seq)
    print(seq_int == [1, 4, 4, 1, 3, 2, 1, 3, 2, 2])  # OK

    seq_int = appendSentinelNumbers(seq_int)
    print(seq_int == [1, 4, 4, 1, 3, 2, 1, 3, 2, 2, 0, 0, 0])  # OK

    pos_table = getPositionTable(seq_int)
    print(pos_table == [1, 4, 7, 10, 2, 5, 8])  # OK

    tuple_list = getTripletListIndexed(seq_int, pos_table)
    print(
        tuple_list
        == [
            ((4, 4, 1), 1),
            ((3, 2, 1), 4),
            ((3, 2, 2), 7),
            ((0, 0, 0), 10),
            ((4, 1, 3), 2),
            ((2, 1, 3), 5),
            ((2, 2, 0), 8),
        ]
    )  # OK

    print(
        getTripletList(tuple_list)
        == [
            (4, 4, 1),
            (3, 2, 1),
            (3, 2, 2),
            (0, 0, 0),
            (4, 1, 3),
            (2, 1, 3),
            (2, 2, 0),
        ]
    )  # OK

    sorted_triplet_list = radixSort(tuple_list)
    print(
        sorted_triplet_list
        == [
            ((0, 0, 0), 10),
            ((2, 1, 3), 5),
            ((2, 2, 0), 8),
            ((3, 2, 1), 4),
            ((3, 2, 2), 7),
            ((4, 4, 1), 1),
            ((4, 1, 3), 2),
        ]
    )

    index_list = getIndexTable(tuple_list)
    print(index_list == [1, 4, 7, 10, 2, 5, 8])  # OK

    index_list_sorted = getIndexTable(sorted_triplet_list)
    print(index_list_sorted == [10, 5, 8, 4, 7, 1, 2])  # OK

    order_table = getOrderTable(sorted_triplet_list)
    print(order_table == [1, 2, 3, 4, 5, 6, 7])  # OK

    print(
        createNextIntList(index_list, index_list_sorted, order_table)
        == [6, 4, 5, 1, 7, 2, 3]
    )  # OK
