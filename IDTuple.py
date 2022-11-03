class IDTuple:
    integer: int
    ID: int

    def __init__(self, integer: int, ID: int) -> None:
        self.integer = integer
        self.ID = ID

    def __str__(self) -> str:
        return f"({self.integer},{self.ID})"

    def isLinked(self, other_IT) -> bool:
        return self.ID == other_IT.ID


if __name__ == "__main__":
    A = IDTuple(56, 0)
    B = IDTuple(12, 0)
    C = IDTuple(776, 1)
    D = IDTuple(59, 2)
    E = IDTuple(93, 3)
    print(A)  # OK
    print(A.isLinked(B) is True)  # OK
    print(A.max(B, C, D, E) == 776)  # OK
