import numpy as np
from Bio import SeqIO
# import cProfile
import bwt
from RESULTS_LIONEL.Import_DC3_Lionel import import_file

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
    print(f"First letter: {firstBase}")

    bottom = 0
    top = lenBwt - 1  # Stay in the string boundaries
    i = len(kmer) - 1  # Stay in the kmer boundaries
    print(f"first e: {bottom},  first f: {top}")

    nbOccur = 1

    while nbOccur > 0 and i >= 0:
        X = kmer[i]

        if X not in bwtDna:
            return False
        else:
            print(f"i: {i}, X: {X}")
            print(f"Rank table of {X}: {rankMat[X]}")

            firstOcc = bwtSort.index(X) - 1
            print(f"First index of {X} in bwtSort: {firstOcc}")

            # Counts all the occurences of X in bwtDna between an empty
            # rank and the following rank to deduce the rank of an
            # actual index e
            preBot = bottom
            rankBot = rankMat[X][preBot]
            print(f"rankE before while: {rankBot}")
            countBot = 0  # number of times the letter X is met in bwtDna
            # We scan bwtDna for X letter until we reach
            while rankBot == -1:
                preBot += 1  # e < f, we increase to stay in the boundaries
                if bwtDna[preBot] == X:
                    countBot += 1
                rankBot = rankMat[X][preBot]
            rankBot -= countBot
            print(f"rankE after while: {rankBot}")

            preTop = top
            rankTop = rankMat[X][preTop]
            print(f"rankF before while: {rankTop}")
            countTop = 0
            while rankTop == -1:
                preTop -= 1  # f > e, we decrease to stay in the boundaries
                if bwtDna[preTop] == X:
                    countTop += 1
                rankTop = rankMat[X][preTop]
            rankTop += countTop
            print(f"rankF after while: {rankTop}")

            # The first character of the BWT has a rank of 1 in the
            # rank matrix but it is the first appearance of this
            # character. To take into account this information, we
            # decrease e of 1 in presence of this character, i.e when
            # there is no character before it
            if bottom == 0 and X == firstBase:
                bottom = firstOcc + rankBot - 1
            else:
                bottom = firstOcc + rankBot
            top = firstOcc + rankTop
            print(f"e: {bottom}, f: {top}")

            nbOccur = top - bottom  # quantity of elements between 2 indexes

            if nbOccur == 0:
                return False
            print(f"Number of pattern: {nbOccur}\n")

            # positions of the kmer in the chromosome
            locs = suffixTab[bottom:top]

            i -= 1

            # True if the entire pattern is crossed
            if i == -1:
                return kmer, np.sort(locs, kind="mergesort")

    return kmer


def verification_pattern(chromo, kmer, locs):
    kmerL = kmer.lower()
    print(
        f"Number of kmer in the chromosome: {chromo.count(kmerL)}")
    print(f"Kmer length: {len(kmerL)}")
    i = 1
    for loc in locs:

        print(f"{i:>3} loc: {loc}")
        print(
            f"{chromo[loc -10:loc]}--{chromo[loc:loc + len(kmer)]}--{chromo[loc + len(kmer) +1:loc + len(kmer) + 11]}")
        i += 1


def compare_ndarray(array1, array2):
    resArray = []
    length = len(array1)
    for i in range(length):
        resArray.append(array1[i] == array2[i])
    return resArray


def mapping(chromo, read):
    pass


if __name__ == "__main__":
    chromo = []
    for record in SeqIO.parse("SEQUENCES/P_fal_genome.fna", format="fasta"):
        chromo.append(str(record.seq))
    print(f"chromosome 1: {chromo[0][:100]}, type: {type(chromo[0])}")

    chromo1 = Chromosome("P_fal_chromosome_1", chromo[0])

    DC3_chrom1 = import_file("DC3_chrom1.txt", "DC3")
    print(
        f"type dc3lionel: {type(DC3_chrom1)} type suffix table {type(chromo1.suffix_table)}")
    print(f"Lengths: {len(DC3_chrom1)} {len(chromo1.suffix_table)}")
    print(
        f"dc3Lionel: {DC3_chrom1[:100]} {DC3_chrom1[len(DC3_chrom1) - 100:]}\nsuffix table: {chromo1.suffix_table[:100]}")
    print(
        f"Comparison array: {compare_ndarray(DC3_chrom1, chromo1.suffix_table)[:100]}")
    # print(f"Suffix table: {s_table}")

    # testing for the first chromosome
    bwtChromo1 = bwt.bwt(chromo[0], chromo1.suffix_table)
    # print(f"bwt: {bwtChromo1[:10000]}, type: {type(bwtChromo1)}")
    print(f"Number of A: {bwtChromo1.count('A')}")

    bwtList = list(bwtChromo1)
    # print(f"Sorted BWT list: {sorted(bwtList)[:100]}")

    # cProfile.run("rankMat = bwt.create_rank_mat(bwtChromo1)")
    rankMat = bwt.create_rank_mat(bwtChromo1)
    # print(f"Rank matrix of A: {rankMat['A'][260000:261000]}")

    pattern = "AAAAAAAAAAAAAAAAAAAAAAAAAA"
    kmer, locs = search_kmer_pos(
        bwtChromo1, rankMat, chromo1.suffix_table, pattern
    )

    print(f"Localisations: {locs}, length DNA: {len(chromo1.DNA)}")
    verification_pattern(chromo1.DNA, pattern, locs)
