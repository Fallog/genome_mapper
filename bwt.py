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
    until now. Because its creation is costful, we decide to count only
    one over 32 indexes, the rank checkpoint. The remaining indexes are
    discarded, filled with -1.

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
        # We use int32 because 32 bits is sufficient to store the
        # values we need
        countings = np.zeros(lenDna, dtype=np.int32)
        # We count only one over 32 indexes to be time efficient
        for i in range(lenDna):
            if i % 32 == 0:
                # :(i + 1) because we want to include the first & last char of bwtDna
                countings[i] = bwtDna[:(i + 1)].count(letter)
            else:
                # row is discarded i.e its rank is not computed
                countings[i] = -1
        rkMat[letter] = countings
    return rkMat


if __name__ == "__main__":
    chromo = []
    for record in SeqIO.parse("SEQUENCES/P_fal_genome.fna", format="fasta"):
        chromo.append(str(record.seq))
    # print(f"chromosome 1: {chromo[0]}")

    s_table = Chromosome("P_fal_chromosome_1", chromo[0]).import_dc3_result()
    print(f"Suffix table: {s_table}")

    bwtChromo1 = bwt(chromo[0], s_table)  # testing for the first chromosome
    print(f"bwt: {bwtChromo1[:100]}, type: {type(bwtChromo1)}")

    bwtList = list(bwtChromo1)
    print(f"Sorted BWT list: {sorted(bwtList)[:100]}")

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
