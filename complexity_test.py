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
# for i in range(30000):
#     seq += random.choice(["A", "C", "G", "T"])
# seq += "$"
# st = y.dc3(seq)
# a = bwt.bwt(seq, st)
# print(a[-1])
# b1 = bwt.create_rank_mat(a)
# b2 = bwt.create_rank_mat2(a)

# c = np.where((b1["A"] == b2["A"]) == False)
# d = np.where((b1["C"] == b2["C"]) == False)
# e = np.where((b1["G"] == b2["G"]) == False)

# print(c)
# print(d)
# print(e)

# print(b1["A"][-1])
# print(b2["A"][-1])

# print(b1["C"][-1])
# print(b2["C"][-1])

# print(b1["G"][-1])
# print(b2["G"][-1])


# cProfile.run("""bwt.create_rank_mat(a)""")
# cProfile.run("""bwt.create_rank_mat2(a)""")
