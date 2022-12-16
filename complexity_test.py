import dc3 as y
import cProfile
import random
import dc3_test
import numpy as np
import tools
import bwt
import mapping


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


def test_rank_mat():
    seq = ""
    for i in range(1000000):
        seq += random.choice(["A", "C", "G", "T"])
    seq += "$"
    st = y.dc3(seq)
    bwtseq = bwt.bwt(seq, st)
    # print(bwtseq)
    print("BWT DONE !")
    return bwtseq


def make_rand_seq(lenght: int):
    seq = ""
    for i in range(lenght - 1):
        seq += random.choice(["A", "C", "G", "T"])
    seq += "$"
    return seq


random.seed(188)
read = make_rand_seq(40)
# print(read)
# res = mapping.cut_read_to_kmer(read, 2)
# print(res)
# cProfile.run("mapping.cut_read_to_kmer(read, 500)")

# a = test_rank_mat()
# cProfile.run("""bwt.create_rank_mat(a)""")
pat = "ACA"
sf = y.dc3(read)
print(sf, "sf")
bwt_dna = bwt.bwt(read, sf)
rank_mat = bwt.create_rank_mat(bwt_dna)

res2 = mapping.search_kmer_pos(bwt_dna, rank_mat, sf, pat)
res_theo = [i for i in range(len(read)) if read.startswith(pat, i)]
print("seq :", read)
print("en théorie :", res_theo)

for i in res_theo:
    print(read[i : i + 2])


print("___________")
print("seq :", read)
print("ton res :", res2)

# for i in res2[1]:
#     print(read[i : i + 2])

res3 = mapping.string_search(bwt_dna, pat, rank_mat, sf)
print("___________")
# print("seq :", read)
# print("".join(sorted(bwt_dna)))
# print(bwt_dna)
# print(rank_mat)
print("mon res :", res3)


for i in res3[0]:
    print(read[i : i + 2])
