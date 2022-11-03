from IDTuple import IDTuple


class IDTupleList:
    size: int
    id_tuple_list = [IDTuple]

    def __init__(self, int_list: list[int]) -> None:
        self.size = len(int_list)
        self.id_tuple_list = [0] * self.size
        for index, integer in enumerate(int_list):
            self.id_tuple_list[index] = IDTuple(integer, index)

    def __str__(self) -> str:
        string = "["
        for ID_tuple in self.id_tuple_list:
            string += str(ID_tuple) + ","
        return string[:-1] + "]"

    def max(self) -> int:
        max_int = 0
        for IT in self.id_tuple_list:
            max_int = max(max_int, IT.integer)
        return max_int

    def count(self, target_int: int) -> int:
        return [self.id_tuple_list[i].integer for i in range(self.size)].count(
            target_int
        )


if __name__ == "__main__":
    int_list = [
        98,
        98,
        99,
        0,
        99,
        99,
        97,
    ]
    id_t_list = IDTupleList(int_list)
    print(id_t_list)  # OK
    print(id_t_list.max() == 99)  # OK
    print(id_t_list.count(99) == 3)  # OK
