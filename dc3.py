# import cProfile
# import random
from radixsort import radix_sort


def convert_dna_to_int(dna_seq: str) -> list[int]:
    """Convert a string of character into a list of integers representing each character.
    Because the strings are DNA sequence, the conversion is made with a dictionary mapping each nucleotide with a number.
    A -> 1 ; C -> 2 ; G -> 3 ; T -> 4
    The chosen number respect the order of this 4 letters in the alphabet.

    Args:
        dna_seq (str): sequence of A, C, G, T

    Returns:
        list[int]: representation of dna_seq with an array of integers
    """
    # DNA dictionary faster then ascii conversion
    dna_dict = {"A": 1, "C": 2, "G": 3, "T": 4}
    int_array = [0] * len(dna_seq)  # performance
    for index, nucleotide in enumerate(dna_seq):
        int_array[index] = dna_dict[nucleotide]
    return int_array


def add_sentinel_numbers(input_list: list[int]) -> list:
    """Returns input_list with three 0 appended at the end of the list.

    Args:
        input_list (list[int]): a list of integer

    Returns:
        list: input_list but concatenated with a list of three 0
    """
    return input_list + [0, 0, 0]


def get_position_table(input_list: list[int], wanted_position: int) -> list[int]:
    """Returns an array containing indexes i of input_list argument that satisfy i % 3 = wanted_position and
    also i + 2 < len(input_list).

    Args:
        input_list (list[int]): array of integers
        wanted_position (int): either 0, 1 or 2. Defines the indexes of input_list that will be put in the returned list

    Returns:
        list[int]: contains indexes i of input_list such that i % 3 = wanted_position.
    """
    size = len(input_list)  # performance
    # Implementation with 2 if loop -> simpler and more straightforward.
    return [i for i in range(size) if i % 3 == wanted_position and i + 2 < size]
    # Possible to do it with range(wanted_position, size, 3) and if i + 2 < size.
    # Both method's performance need to be tested


def get_triplet_index(
    input_list: list[int], position_table: list[int]
) -> list[tuple[tuple[int, int, int], int]]:
    """Returns a list of tuple made of a triplet plus its index in the position table.
    The triplets are built using the position table and an array of integers.
    For each element p in the position table, we form a tuple looking like this :
    ((input_list[p], input_list[p + 1], input_list[p + 2]), p)

    Args:
        input_list (list[int]): an array of integers
        position_table (list[int]): array obtained by the get_position_table function

    Returns:
        list[tuple[tuple[int, int, int], int]]: list of tuple where a tuple consists in a triplet and
    an integer representing the index from which the triplet is built
    """
    pos_size = len(position_table)  # performance
    triplet_list = [0] * pos_size  # performance
    for i in range(pos_size):
        position = position_table[i]  # performance
        triplet_list[i] = (
            (
                input_list[position],
                input_list[position + 1],
                input_list[position + 2],
            ),
            position,
        )
    return triplet_list


def get_triplet(
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
    # Intuitive and fast to implement
    return [index_triplet_list[i][0] for i in range(len(index_triplet_list))]


def get_indexes(triplet_list: list[tuple[tuple[int, int, int], int]]) -> list[int]:
    """Returns an array of integers built with the int inside the tuple of triplet_list.
    triplet_list is an array of ((int, int, int), int) and here it extracts the last int and put it in the returned list.

    Args:
        triplet_list (list[tuple[tuple[int, int, int], int]]): a tuple consists in a triplet and an integer representing
        the index from which the triplet was built in the first place

    Returns:
        list[int]: indexes of the triplet in an sorted or unsorted triplet_list
    """
    return [triplet_list[i][1] for i in range(len(triplet_list))]


def get_orders(sorted_tuple_list: list[tuple[tuple[int, int, int], int]]) -> list[int]:
    """Returns the sorted position of each triplet in sorted_tuple_list.
    If two triplet are identical, their position is set equal.

    Args:
        sorted_tuple_list (list[tuple[tuple[int, int, int], int]]): a tuple consists in a triplet and an integer
        representing the index from which the triplet was built in the first place

    Returns:
        list[int]: sorted positions of the triplets
    """
    size = len(tuple_list)  # performance
    # element 0 have order 1 and the for loop start to 1
    order_table = [1] * size

    for i in range(1, size):
        parsed_tuple = sorted_tuple_list[i][0]  # selecting the tuple
        previous_tuple = sorted_tuple_list[i - 1][0]
        if parsed_tuple == previous_tuple:
            order_table[i] = order_table[i - 1]
        else:
            order_table[i] = order_table[i - 1] + 1
    return order_table


def create_next_list(
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
    size = len(index_table)  # performance
    next_int_list = [0] * size  # performance
    for i in range(size):
        next_int_list[i] = order_table[index_table_sort.index(index_table[i])]
    return next_int_list


def get_pairs_index(
    int_sequence: list[int], position_table: list[int], index_table: list[int]
) -> list[tuple[tuple[int, int], int]]:
    size = len(position_table)  # performance
    pairs_table = [0] * size  # performance
    for i in range(size):
        position = position_table[i]
        pairs_table[i] = (
            (
                int_sequence[position],
                index_table[int_sequence[position + 1]
                            - 1] if i != size - 1 else 1,
            ),
            position,
        )
    return pairs_table


def get_smallest_index(int_seq: list[int], index12: int, index0: int, nb_loop: int = 0) -> int:
    """_summary_

    Args:
        int_seq (list[int]): _description_
        index12 (int): _description_
        index0 (int): _description_
        nb_loop (int, optional): _description_. Defaults to 0.

    Returns:
        int: _description_
    """
    if int_seq[index12] == int_seq[index0]:
        return get_smallest_index(int_seq, index12 + 1, index0 + 1, nb_loop + 1)
    elif int_seq[index12] > int_seq[index0]:
        return index0 - nb_loop
    else:
        return index12 - nb_loop


def merge_index_tables(
    int_seq: list[int], index12_table: list[int], index0_table: list[int]
) -> list[int]:
    """_summary_

    Args:
        int_seq (list[int]): _description_
        index12_table (list[int]): _description_
        index0_table (list[int]): _description_

    Returns:
        list[int]: _description_
    """
    size12, size0 = len(index12_table), len(index0_table)
    index_table_merged = []  # because table is not filled in one loop
    i, j = 0, 0  # one increment for each table i for
    # The while loop stops when one of the increment went throug all its list
    while i < size12 and j < size0:
        # Next 2 lines are used for better readability
        index12 = index12_table[i]
        index0 = index0_table[j]
        smallest_index = get_smallest_index(int_seq, index12, index0)
        index_table_merged.append(smallest_index)
        if smallest_index in index12_table:
            i += 1
        else:
            j += 1
    # When the while loop stops, we add the remaining list to index_table_merged
    return index_table_merged + index12_table[i:] + index0_table[j:]


def dc3(int_sequence: list[int]):
    pass


if __name__ == "__main__":
    test_seq = "ATTAGCAGCC"
    print("Length of the sequence: ", len(test_seq))
    seq_int = convert_dna_to_int(test_seq)
    print(seq_int == [1, 4, 4, 1, 3, 2, 1, 3, 2, 2])  # OK

    seq_int = add_sentinel_numbers(seq_int)
    print(seq_int == [1, 4, 4, 1, 3, 2, 1, 3, 2, 2, 0, 0, 0])  # OK

    pos_table_1 = get_position_table(seq_int, wanted_position=1)
    print(pos_table_1 == [1, 4, 7, 10])  # OK
    pos_table_2 = get_position_table(seq_int, wanted_position=2)
    pos_table = pos_table_1 + pos_table_2
    print(pos_table == [1, 4, 7, 10, 2, 5, 8])  # OK

    tuple_list = get_triplet_index(seq_int, pos_table)
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
        get_triplet(tuple_list)
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

    sorted_triplet_list = radix_sort(tuple_list)
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

    index_list = get_indexes(tuple_list)
    print(index_list == [1, 4, 7, 10, 2, 5, 8])  # OK

    index_list_sorted = get_indexes(sorted_triplet_list)
    print(index_list_sorted == [10, 5, 8, 4, 7, 1, 2])  # OK

    order_table = get_orders(sorted_triplet_list)
    print(order_table == [1, 2, 3, 4, 5, 6, 7])  # OK

    next_int_list = add_sentinel_numbers(
        create_next_list(index_list, index_list_sorted, order_table)
    )
    print(next_int_list == [6, 4, 5, 1, 7, 2, 3, 0, 0, 0])  # OK

    next_pos_table_0 = get_position_table(next_int_list, wanted_position=0)
    print(next_pos_table_0 == [0, 3, 6])  # OK

    pairs_table = get_pairs_index(
        int_sequence=[3, 3, 4, 1, 4, 5, 2, 0, 0, 0],
        position_table=[0, 3, 6],
        index_table=[7, 1, 2, 4, 5],
    )
    print(pairs_table == [((3, 2), 0), ((1, 4), 3), ((2, 1), 6)])  # OK

    sort_pairs_table = radix_sort(pairs_table)
    index0 = get_indexes(sort_pairs_table)
    print(sort_pairs_table, index0)

    print(get_smallest_index([3, 3, 4, 1, 4, 5, 2, 0, 0, 0], 7, 3) == 7)  # OK
    print(get_smallest_index([3, 3, 4, 1, 4, 5, 2, 0, 0, 0], 1, 3) == 3)  # OK
    print(get_smallest_index([3, 3, 4, 1, 4, 5, 2, 0, 0, 0], 1, 0) == 0)  # OK
    print(get_smallest_index(
        [97, 98, 99, 97, 98, 99, 97, 99, 97, 98, 0, 0, 0], 8, 0) == 8)

    print(merge_index_tables([3, 3, 4, 1, 4, 5, 2, 0, 0, 0], [7, 1, 2, 4, 5], [3, 6, 0])
          == [7, 3, 6, 0, 1, 2, 4, 5])  # OK
    print(merge_index_tables([97, 98, 99, 97, 98, 99, 97, 99, 97, 98, 0, 0, 0],
                             [10, 8, 1, 4, 7, 2, 5], [0, 3, 6, 9]) == [10, 8, 0, 3, 6, 9, 1, 4, 7, 2, 5])  # OK
