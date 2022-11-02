# import numpy   ## A convertir en numpy
from array import array
import numpy as np
from quickSort import quickSort


def strToAscii(text: str) -> array:  # ya un numpy loadtxt qui existe
    res_array = np.frombuffer(text.encode(), np.int8)
    return res_array


def asciiToStr(ascii_list: array) -> str:
    res_txt = ascii_list.tobytes().decode("ascii")
    return res_txt


def makeP12(t_list: array) -> array:
    lenght = len(t_list)
    p1 = np.arange(1, lenght - 2, 3)
    p2 = np.arange(2, lenght - 2, 3)
    return np.concatenate((p1, p2))


def makeTriplet(t_list: list[int], p12_list: list[int]) -> list[list[int]]:
    r12 = np.zeros(len(p12_list), dtype=[("1", int), ("2", int), ("3", int)])
    for i, ind in enumerate(p12_list):
        r12[i] = (t_list[ind], t_list[ind + 1], t_list[ind + 2])
    return r12


def sortTriplet(r12: list[list[int]]) -> list[list[int]]:
    return quickSort(r12)


def occur_triplet_list(r12_sorted: list[list[int]]) -> list[int]:
    trash_count = []
    counter = 0
    for i in r12_sorted:
        if i in trash_count:
            pass


def step1(init_seq):
    T = np.concatenate((init_seq, np.zeros(3)))
    p12 = makeP12(T)
    print(p12)
    r12 = makeTriplet(T, p12)
    print(r12)
    tp = sortTriplet(r12)
    print(tp)
    print(r12)

    if tp.shape == np.unique(tp).shape:
        return None  # Jcp pas encore ce qu'on doit faire apres
    else:
        step1(tp)


if __name__ == "__main__":
    s = "abcabcacab"

    T = strToAscii(s)

    step1(T)


# Important pas stocker tout les info, en théorie y'a juste order a stocker. (P0 a voir ptete besoin a la fin).
# Recoder le "radix" short algo
