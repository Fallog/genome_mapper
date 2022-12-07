import numpy as np
from Chromosome import Chromosome
from Bio import SeqIO


def bwt(dna, suffixTable, end_of_string="$"):
    """Returns the Burrows Wheeler Transform of the dna argument.
    The BWT is computed efficiently using the suffix array of dna.

    Args:
        dna (str): string made only with A, T, C and G character
        suffix_table (ndarray): contains all the suffixes of dna
        end_of_string (char): appended character specifying the end of
            the string

    Return:
        str: BWT transformation of dna
    """
    bwtStr = ""
    dna += end_of_string

    for i in range(len(suffixTable)):
        # Chromosomes have sometimes bases in lowercase
        bwtStr += dna[suffixTable[i] - 1].capitalize()

    return bwtStr


def create_rank_mat(bwtDna):
    """Returns a dict storing the rank tables of each character in
    the dnaSeq argument.
    A rank table is an array containing, for each character of bwtDna
    and at each position, the number of occurences of the character
    until now.

    Args:
        dnaSeq (str): a string made only with A, T, C, G and a $

    Returns:
        dict: contains 5 keys: A, T, C, G and $, each key is mapped
            with the rank table in the dnaSeq of the associated
            character
    """
    alphabet = ['A', 'T', 'C', 'G', '$']
    rkMat = {}
    lenDna = len(bwtDna)
    for letter in alphabet:
        # We use uint32 because 32 bits is sufficient to store the
        # unsigned values we need
        countings = np.zeros(lenDna, dtype=np.uint32)
        # We count only one over 32 rows to be time efficient
        for i in range(lenDna):
            if i % 32 == 0:
                # :(i + 1) because we want to include the first & last char of string
                countings[i] = bwtDna[:(i + 1)].count(letter)
            else:
                countings[i] = 0
        rkMat[letter] = countings
    return rkMat


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
        int: number of times kmer can be read on the DNA sequence
        ndarray: every localisation where kmer can be read
    """
    # Avoid modifying bwtGen while sorting
    bwtSort = sorted(list(bwtDna))
    lenBwt = len(bwtDna)
    # print(f"Sorted BWT: {bwtSort}")

    e = 0
    f = lenBwt - 1  # Stay in the string boundaries
    i = len(kmer) - 1  # Stay in the kmer boundaries
    print(f"first e: {e},  f: {f}")

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
            countE = 0
            while rankE == 0:
                newE += 1
                if bwtDna[newE] == X:
                    countE += 1
                rankE = rankMat[X][newE]
            rankE -= countE
            print(f"rankE: {rankE}")

            newF = f
            rankF = rankMat[X][newF]
            countF = 0
            while rankF == 0:
                newF -= 1
                if bwtDna[newF] == X:
                    countF += 1
                rankF = rankMat[X][newF]
            rankF -= countF
            print(f"rankF: {rankF}")

            # The first character of the BWT has a rank of 1 in the
            # rank matrix but it is the first appearance of this
            # character. To take into account this information, we
            # decrease e of 1 in presence of this character, i.e when
            # there is no character before it
            if len(bwtDna[:(e + firstOcc)]) == 0:
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
                return kmer, nbOccur, locs

    return kmer


if __name__ == "__main__":
    chromo = []
    for record in SeqIO.parse("SEQUENCES/P_fal_genome.fna", format="fasta"):
        chromo.append(str(record.seq))
    # print(f"chromosome 1: {chromo[0]}")

    s_table = Chromosome("P_fal_chromosome_1").import_dc3_result()
    # print(f"Suffix table: {s_table}")

    bwtChromo1 = bwt(chromo[0], s_table)  # testing for the first chromosome
    # print(f"bwt: {bwtChromo1[:100]}, type: {type(bwtChromo1)}")

    bwtList = list(bwtChromo1)
    # print(f"Sorted BWT list: {sorted(bwtList)[:100]}")

    rankMat = create_rank_mat(bwtChromo1)
    print(f"Rank matrix of A: {rankMat['A'][260000:261000]}")

    # disp = "[ "
    # for count in rankMat['A']:
    #     disp += f" {count}"
    # disp += "]"
    # print(disp)

    # bwtT = bwt(T)
    # print(f"BWT from suffix table: {bwtT}")

    # bwtDC3 = bwt(T)
    # print(f"Is the BWT DC3 the same ? {bwtDC3 == bwtT}")
    # print(f"BWT DC3: {bwtDC3}")

    print(search_kmer_pos(bwtChromo1, rankMat, s_table, "ATA"))
