import dc3 as y
import cProfile
import random
import dc3_test
import numpy as np
import tools
import bwt


def test_dc3():
    # seq = "AGAGAGGAGCGCGCGAGAGCGCGAGAG" * 500
    # seq += "$"
    # random.seed(12)
    seq = ""
    for i in range(600000):
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
    # print(c[0].size)
    # print(a == b)
    print(c)


seq = ""
for i in range(15):
    seq += random.choice(["A", "C", "G", "T"])
seq += "$"
import mapping

print(seq)
print(len(seq))
print(mapping.cut_read_to_kmer(seq, 10))

# cProfile.run("""mapping.cut_read_to_kmer("seq", 10)""")
