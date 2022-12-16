import numpy as np
from Bio import SeqIO
from tqdm import tqdm
from Chromosome import Chromosome


def verification_pattern(dna_seq: str, pattern: str, locs: np.ndarray):
    """Verifies the founded localisation of the pattern argument over
    the dna_seq argument.

    Args:
        dna_seq (str): big string made only with A, T, C and G character
        pattern (str): shorter string made with the same letters
        locs (ndarray): every founded localisation of pattern over
            dna_seq
    """
    if locs[0] == -1:
        print("No localisation is found")
    pattern_L = pattern.lower()  #
    print(f"Kmer length: {len(pattern_L)}")
    for i, loc in enumerate(locs):
        print(f"{i:>3} loc: {loc}")
        print(
            f"{dna_seq[loc -10:loc]}--{dna_seq[loc:loc + len(pattern)]}--{dna_seq[loc + len(pattern) +1:loc + len(pattern) + 11]}"
        )


def get_first_occ(rank_mat: dict[np.ndarray], base: str) -> int:
    """Returns the index of the first appearance of the base argument
    in the sorted version of the bwt_dna.

    Args:
        rank_table (dict[np.ndarray]): each ndarray stores the ranks of
            each nucleotide
        base (str): nucleotide for which we want its first appearance

    Returns:
        int: index of the first appearance of base argument
    """
    if base == "A":
        return 0
    elif base == "C":
        return rank_mat["A"][-1]
    elif base == "G":
        return rank_mat["A"][-1] + rank_mat["C"][-1]
    else:
        return rank_mat["A"][-1] + rank_mat["C"][-1] + rank_mat["G"][-1]


def string_search(bwt: str, read: str, rank_mat: dict[np.ndarray], suff_t: np.ndarray) -> np.ndarray:
    """Returns all the localisation of read argument over the DNA
    sequence that gave the bwt argument. The rank_mat is used to fasten
    the search and the suff_t argument gives the actual localisation.

    Args:
        bwt (str): big string made only with A, T, C and G character
        read (str): shorter string made only with A, T, C and G character
        rank_mat (dict[np.ndarray]): each ndarray stores the ranks of
            each nucleotide
        suff_t (np.ndarray): contains all the suffixes of the DNA
            sequence used to build bwt

    Returns:
        np.ndarray: array of the localisations of read over the DNA 
    """
    top = 0
    bottom = len(bwt) - 1
    while top <= bottom and read != "":
        base = read[-1]
        read = read[: len(read) - 1]
        top = get_first_occ(rank_mat, base) + rank_mat[base][top]
        if top < bottom:
            bottom = (
                get_first_occ(rank_mat, base) + rank_mat[base][bottom + 1] - 1
            )  # Because it will count 2 times the first
        else:
            bottom = get_first_occ(rank_mat, base) + rank_mat[base][bottom]
    if read == "":
        return suff_t[top: bottom + 1]
    else:
        return np.array([-1])


def cut_read_to_kmer(read: str, patt_len: int) -> list[str]:
    """Divide the given read argument into a list of k-mer,
    smaller strings of size the k_len argument.

    Args:
        read (str): str made of the characters 'A', 'T', 'C' and 'G'
        patt_len (int): size of the k-mer to be created from the read
            sequence

    Returns:
        list[str]: list of all the k-mer created from the read
    """
    return [read[i: i + patt_len] for i in range(0, len(read), patt_len)]


def get_read_quality(first_loc: np.ndarray, loc_list: np.ndarray, patt_len: str):
    """Returns the "quality" of a read, i.e the number of localisation
    it is possible to retrieve between the 2 extreme localisations
    (first and last indexes of loc_list)


    Args:
        first_loc (np.ndarray): array of every localisations of the
            first pattern
        loc_list (np.ndarray): array of every localisation for each
            pattern of a given read
        patt_len (np.ndarray): length of the pattern a read is divided
            into

    Returns:
        int: number of mismatched patterns in between first and last
            patterns in the patterns list 
    """
    targ_loc = [first_loc + (i * patt_len) for i in range(1, patt_len)]
    # print(f"Locs we are looking for: {targetLoca}")
    # Counter giving the number of unplaced kmer between the extreme reads
    nb_mismatch = 0
    index_patt = 1  # Index specifying on which kmer we are
    index_parsed = 0  # Index to parse the localisations of inside locaList
    index_targ = 0  # Index to parse targetLoca
    while index_patt < patt_len - 1:
        locaParsed = loc_list[index_patt]
        # print(f"Locations parsed: {locaParsed} len: {len(locaParsed)}")
        # If it reaches the end of actual localisation list
        # before returning anything go to next kmer but count one
        # unplaced kmer
        if index_parsed == len(locaParsed):
            # print(f"indParsed: {indParsed}")
            index_patt += 1
            index_targ += 1
            index_parsed = 0
            nb_mismatch += 1
        else:
            # If the positions match then go to next kmer and next target
            # localisation
            if locaParsed[index_parsed] == targ_loc[index_targ]:
                # print("Found !")
                index_patt += 1
                index_targ += 1
                index_parsed = 0
            else:
                index_parsed += 1
    return nb_mismatch


def link_kmer_fast(patt_list, loc_lis):
    """Returns the reconstructed read and its localisation on a
    chromosome.
    Linkage of the 2 extreme kmers of kmerList (first and last) consists
    in finding if the localisation of the last is equal to the
    localisation of the first + N-1 times the length of a kmer, (N being
    the number of kmer in kmerList).

    Args:
        kmerlist (list[str]): list of every kmer of a given read, output
            of cut_read_to_kmer function
        locaList (ndarray): array of every localisation for each kmer
            in kmerList, output of search_kmer_pos

    Returns:
        str: read built from all the kmer in kmerList
        ndarray: localisation of the read on a given chromosome,
            np.empty(1) if no localisation is found
    """
    patt_nb = len(patt_list)
    patt_len = len(patt_list[0])
    first_loc = loc_lis[0]
    last_loc = loc_lis[-1]
    first_index = 0
    last_index = 0
    read = "".join(patt_list)
    # If first or last kmer have 0 localisation, return 0 localisation
    # for the read
    if len(last_loc) == 1 and last_loc[0] < 1:
        # print("No kmer localisation")
        return read, np.empty(1)
    if len(first_loc) == 1 and first_loc[0]:
        # print("No kmer localisation")
        return read, np.empty(1)
    locaRead = []
    while first_index != len(first_loc):
        if last_index == len(last_loc):
            # print("Next first localisation")
            first_index += 1
            last_index = 0
        else:
            expLoca = first_loc[first_index] + patt_len * (patt_nb - 1)
            if expLoca >= last_loc[last_index]:
                if expLoca == last_loc[last_index]:
                    # print("Equal")
                    locaRead.append(first_loc[first_index])
                    first_index += 1
                    last_index = 0
                else:
                    # print("Superior")
                    last_index += 1
            else:
                # wprint("Inferior")
                first_index += 1
                last_index = 0
    return read, np.array(locaRead)


def mapping(chromo, read):
    pass


if __name__ == "__main__":
    chromo = []
    with tqdm(total=15, desc="Chromosomes importation") as pbar:
        for record in SeqIO.parse("SEQUENCES/P_fal_genome.fna", format="fasta"):
            chromo.append(str(record.seq))
            pbar.update(1)
    chromo1 = Chromosome("P_fal_chromosome_1", chromo[0], 1)

    # AAACCCTGAA$
    # seq = "ATAATA$"
    # suffixT = dc3(seq=seq)
    # print(f"suffix array: {suffixT} ST[1:4]: {suffixT[1:4]}")
    # bwtSeq = bwt.bwt(seq, suffixT)
    # print(f"bwt: {bwtSeq}")
    # rkMtSeq = bwt.create_rank_mat(bwtSeq)
    # print(f"Sorted bwt: {sorted(list(bwtSeq))}")
    # locs = string_search(bwtSeq, "AT", rkMtSeq, suffixT)
    # print(locs, "tata")
    # verification_pattern(seq, "AT", locs)

    # reads = []
    # with tqdm(total=1500000, desc="Reads importation") as pbar:
    #     for record in tqdm(SeqIO.parse("SEQUENCES/P_fal_reads.fq", format="fastq")):
    #         reads.append((record.seq, record.id))
    #         pbar.update(1)
    # # First read of mapping_P_fal.sam file
    # readTest = str(reads[100000 - 21515][0])
    # print(f"First read: {readTest} type: {type(readTest)}")
    readTest = "AAACCCTGAACCCTAAACCCTGAACCCTAAACCCTAAACCCTGAACCCTAAACCCTAAACCCTGAACCCTAAACCCTGAAACCTAAAACCTGAACCCTAA"
    kmerFstRead = cut_read_to_kmer(readTest, 10)
    print(f"Kmers: {kmerFstRead}")

    # testing for the first chromosome
    bwtChromo1 = chromo1.bwt
    # print(f"bwt: {bwtChromo1[:10000]}, type: {type(bwtChromo1)}")

    rankMat = chromo1.rank_mat
    # print(f"First occ of C {get_first_occ(rankMat, 'C')}")  # OK
    # print(f"Rank matrix of C: {rankMat['C'][:100]}")

    locs = string_search(
        bwtChromo1, kmerFstRead[0], rankMat, chromo1.suffix_table)
    # for kmer in kmerFstRead:
    #     locs.append(search_kmer_pos(bwtChromo1, rankMat,
    #                 chromo1.suffix_table, kmer)[1])
    # print(f"Kmer localisation: {locs}")
    # verification_pattern(chromo1.DNA, kmerFstRead[0], locs)
    # recoReadFst = link_kmer_fast(kmerFstRead, locs)
    # print(f"Reconstructed read: {recoReadFst}")
    # qltyPos1 = get_read_quality(recoReadFst[1][0], locs, 10)  # 1
    # qltyPos2 = get_read_quality(recoReadFst[1][1], locs, 10)  # 5
    # print(f"Type output link_kmer: {type(recoReadFst[1][0])}")
    # print(f"Read qlty1: {qltyPos1} 2: {qltyPos2}")

    bwtList = list(bwtChromo1)
    # print(f"Sorted BWT list: {sorted(bwtList)[:100]}")

    # cProfile.run("rankMat = bwt.create_rank_mat(bwtChromo1)")
    rankMat = chromo1.rank_mat.item()
    # print(f"Rank matrix of A: {rankMat['A'][260000:261000]}")

    locs = []
    for kmer in kmerFstRead:
        locs.append(search_kmer_pos(bwtChromo1, rankMat,
                    chromo1.suffix_table, kmer)[1])
    print(f"Kmers localisation: {locs}")
    verification_pattern(chromo1.DNA, kmerFstRead[0], locs[0])
    recoReadFst = link_kmer_fast(kmerFstRead, locs)
    print(f"Reconstructed read: {recoReadFst}")
    qltyPos1 = get_read_quality(recoReadFst[1][0], locs, 10)  # 1
    qltyPos2 = get_read_quality(recoReadFst[1][1], locs, 10)  # 5
    print(f"Type output link_kmer: {type(recoReadFst[1][0])}")
    print(f"Read qlty1: {qltyPos1} 2: {qltyPos2}")
    # print(f"Reconstructed read: {recoRead} Lenght read: {len(recoRead[0])}")
    # print(f"Actual read: {readTest:>7}")
