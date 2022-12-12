import y_dc3 as y
import dc3 as s
import cProfile
import random
import dc3_test
import numpy as np
from Sais import Sais
import tools

# seq = "AGAGAGGAGCGCGCGAGAGCGCGAGAG" * 500
# seq += "$"
# random.seed(12)
seq = ""
for i in range(4000000):
    seq += random.choice(["A", "C", "G", "T"])
seq += "$"
# seq = "ACGTGCCTAGCCTACCGTGAAACCGGCGGCGGGAAGCGCGAAAGAGAGAGCGGAGAGAGGCGCGCC$"
# print(seq)
# print(y.dc3(seq))

cProfile.run("""y.dc3(seq, sorting_algorithm = "stable")""")
# print(y.dc3(seq, Tru2
a = y.dc3(seq, "mergesort")

## SAIS TEST : PAS OPTI FINALLEMENT
# test = Sais()
# d = tools.strToBase(seq) - 1
# a = test.build_suffix_array(d)
# cProfile.run("""test.build_suffix_array(d)""")

cProfile.run("""dc3_test.skew(seq)""")
b = dc3_test.skew(seq)

c = np.where((a == b) == False)

print(c)
# print(c[0].size)
# print(a == b)


## 44 sec : 600 000
## 5 minutes : 3M3
