# import numpy   ## A convertir en numpy
from array import array
import numpy as np
from radixSort import radixSort


def strToAscii(text: str) -> array:  # ya un numpy loadtxt qui existe
    res_array = np.fromstring(text, np.int8)
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
    r12 = np.empty((len(p12_list), 3))
    dtype = [("1", int), ("2", int), ("3", int)]
    print(r12)
    for i, ind in enumerate(p12_list):
        r12[i] = np.array((t_list[ind], t_list[ind + 1], t_list[ind + 2]), dtype=dtype)
    return r12


def sortTriplet(r12: list[list[int]]) -> list[list[int]]:
    pass


def occur_triplet_list(r12_sorted: list[list[int]]) -> list[int]:
    trash_count = []
    counter = 0
    for i in r12_sorted:
        if i in trash_count:
            pass


# Idée --> pour sort les triplers : faire en sort qu'un triplé = un nombre a 3*x chiffre.


s = "abcabcacab"

T = strToAscii(s)
T = np.concatenate((T, np.zeros(3)))
print(T)

p12 = makeP12(T)
print(p12)

r12 = makeTriplet(T, p12)
print(r12)

# Important pas stocker tout les info, en théorie y'a juste order a stocker. (P0 a voir ptete besoin a la fin).
# Recoder le "radix" short algo
