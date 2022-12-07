import numpy as np
from Bio import SeqIO
# import cProfile
import bwt

from Chromosome import Chromosome


def search_kmer_pos(bwtDna, rankMat, suffixTab, kmer):
    """Search a kmer in a Burrows Wheeler tranformed chromosome.
    Returns the kmer and it's position(s) in the genome. If the kmer is
    absent, it returns the kmer alone.

    Args:
        bwtDna (str): the BWT of a string made only with A, T, C, G
            and a $
        rankMat (dict): contains 5 keys: A, T, C, G and $, each key is
            mapped with the rank table in the dnaSeq of the associated
            character
        suffixTab (ndarray): sorted array of all the suffixes
            used to built bwtDna
        kmer (str): nucleotides pattern to search in the genome

    Return:
        str: kmer argument
        ndarray: every localisation where kmer can be read
    """
    bwtSort = sorted(list(bwtDna))
    lenBwt = len(bwtDna)
    # print(f"Sorted BWT: {bwtSort}")

    firstBase = bwtDna[0]

    e = 0
    f = lenBwt - 1  # Stay in the string boundaries
    i = len(kmer) - 1  # Stay in the kmer boundaries
    print(f"first e: {e},  first f: {f}")

    nbOccur = 1

    while nbOccur > 0 and i >= 0:
        X = kmer[i]

        if X not in bwtDna:
            return False
        else:
            print(f"i: {i}, X: {X}")
            print(f"Rank table of {X}: {rankMat[X]}")

            firstOcc = bwtSort.index(X) - 1  # $ isn't taken into account
            print(f"First index of {X} in bwtSort: {firstOcc}")

            # Counts all the occurences of X in bwtDna between an empty
            # rank and the following rank to deduce the rank of an
            # actual index e
            newE = e
            rankE = rankMat[X][newE]
            print(f"rankE before while: {rankE}")
            countE = 0  # number of times the letter X is met in bwtDna
            # We scan bwtDna for X letter until we reach
            while rankE == -1:
                newE += 1  # e < f, we increase to stay in the boundaries
                if bwtDna[newE] == X:
                    countE += 1
                rankE = rankMat[X][newE]
            rankE -= countE
            print(f"rankE after while: {rankE}")

            newF = f
            rankF = rankMat[X][newF]
            print(f"rankF before while: {rankF}")
            countF = 0
            while rankF == -1:
                newF -= 1  # f > e, we decrease to stay in the boundaries
                if bwtDna[newF] == X:
                    countF += 1
                rankF = rankMat[X][newF]
            rankF += countF
            print(f"rankF after while: {rankF}")

            # The first character of the BWT has a rank of 1 in the
            # rank matrix but it is the first appearance of this
            # character. To take into account this information, we
            # decrease e of 1 in presence of this character, i.e when
            # there is no character before it
            if e == 0 and X == firstBase:
                e = firstOcc + rankE - 1
            else:
                e = firstOcc + rankE
            f = firstOcc + rankF
            print(f"e: {e}, f: {f}")

            nbOccur = f - e  # quantity of elements between 2 indexes

            if nbOccur == 0:
                return False
            print(f"Number of pattern: {nbOccur}\n")

            locs = suffixTab[e:f]  # positions of the kmer in the chromosome

            i -= 1

            # True if the entire pattern is crossed
            if i == -1:
                return kmer, np.sort(locs, kind="mergesort")

    return kmer


def mapping(chromo, read):
    pass


if __name__ == "__main__":
    chromo = []
    for record in SeqIO.parse("SEQUENCES/P_fal_genome.fna", format="fasta"):
        chromo.append(str(record.seq))
    print(f"chromosome 1: {chromo[0][:100]}, type: {type(chromo[0])}")

    with open("chromosome1.txt", mode='w') as file:
        file.write(chromo[0])

    s_table = Chromosome("P_fal_chromosome_1", chromo[0]).import_dc3_result()
    # print(f"Suffix table: {s_table}")

    # testing for the first chromosome
    bwtChromo1 = bwt.bwt(chromo[0], s_table)
    # print(f"bwt: {bwtChromo1[:10000]}, type: {type(bwtChromo1)}")
    print(f"Number of A: {bwtChromo1.count('A')}")

    bwtList = list(bwtChromo1)
    # print(f"Sorted BWT list: {sorted(bwtList)[:100]}")

    # cProfile.run("rankMat = bwt.create_rank_mat(bwtChromo1)")
    rankMat = bwt.create_rank_mat(bwtChromo1)
    # print(f"Rank matrix of A: {rankMat['A'][260000:261000]}")

    print(search_kmer_pos(bwtChromo1, rankMat,
          s_table, "AAAAAAAAAAAAAAAAAAAAAAAAAA"))
