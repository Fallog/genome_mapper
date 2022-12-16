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
    readLen = len(read)  # Performance
    return [read[i: i + patt_len] for i in range(0, readLen, patt_len)]


def link_kmer(kmerList, locaList):
    """_summary_

    Args:
        kmerList (list[str]): _description_
        locaList (list[ndarray]): _description_

    Returns:
        _type_: _description_
    """
    indFirst = 0  # index to parse the localisation of the first kmer
    indLoca = 0  # index to parse the localisation of the following kmers
    indKmer = 1  # index specifying on which following kmer we are
    locaFirst = locaList[0]  # localisations list of the first kmer
    locaNext = locaList[indKmer]  # localisations list of the following kmer
    read = kmerList[0]
    # bestRead = read
    # bestLoca = locaFirst[indFirst]
    readLen = len(kmerList)  # number of kmer in a read
    kmerLen = len(kmerList[0])  # kmer have all the same size
    linkedKmer = 0
    while indKmer != readLen:
        print(f"Length locaFirst: {len(locaFirst)} locaNext: {len(locaNext)}")
        # If it reaches the last localisation of the first kmer
        # it
        if indFirst == len(locaFirst):
            return read, locaFirst[indFirst] + 1
        else:
            # If the next kmer don't have any localisation
            if len(locaNext) == 1 and locaNext[indLoca] < 1:
                read += "-" * len(kmerList[0])
                indKmer += 1
            else:
                # Si on arrive au bout des localisations, on passe à la
                # localisation suivante du premier kmer en retournant à
                # la première localisation du kmer suivant (le 1 donc)
                if indLoca == len(locaNext):
                    if linkedKmer == 0:
                        print("No more kmer")
                        indFirst += 1
                        indKmer = 1
                        indLoca = 0
                        read = kmerList[0]
                    else:
                        indLoca += 1
                        read += "-" * len(kmerList[0])
                    linkedKmer = 0
                else:
                    print(f"Loca first kmer: {locaFirst[indFirst]}")
                    print(f"Loca parsed kmer: {locaNext[indLoca]}")

                    # Si les localisations se suivent, passer au kmer suivant
                    # tout en revenant à sa première localisation
                    print(
                        f"indFirst: {indFirst} indKmer: {indKmer} indLoca: {indLoca}")
                    if locaFirst[indFirst] + kmerLen * indKmer == locaNext[indLoca]:
                        print("Equal")
                        read += kmerList[indKmer]
                        linkedKmer += 1
                        indKmer += 1
                        indLoca = 0

                    else:  # sinon regarder la localisation suivante du kmer actuel
                        if locaFirst[indFirst] + kmerLen * indKmer > locaNext[indLoca]:
                            print("Superior")
                            indLoca += 1
                        else:
                            print("Inferior")
                            indFirst += 1
                            indKmer = 1
                            indLoca = 0
                            read = kmerList[0]

        # Update if needed the returned read
        # if len(read) > len(bestRead):
        #     bestRead = read
        #     bestLoca = locaFirst[indFirst]
        if len(read) == kmerLen * len(kmerList):
            # First base of the genome is set to base 1 and not 0
            return read, locaFirst[indFirst] + 1
        locaNext = locaList[indKmer]


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
    targetLoca = [first_loc + (i * patt_len) for i in range(1, patt_len)]
    # print(f"Locs we are looking for: {targetLoca}")
    # Counter giving the number of unplaced kmer between the extreme reads
    mismatchKmer = 0
    indKmer = 1  # Index specifying on which kmer we are
    indParsed = 0  # Index to parse the localisations of inside locaList
    indTarg = 0  # Index to parse targetLoca
    while indKmer < patt_len - 1:
        locaParsed = loc_list[indKmer]
        # print(f"Locations parsed: {locaParsed} len: {len(locaParsed)}")
        # If it reaches the end of actual localisation list
        # before returning anything go to next kmer but count one
        # unplaced kmer
        if indParsed == len(locaParsed):
            # print(f"indParsed: {indParsed}")
            indKmer += 1
            indTarg += 1
            indParsed = 0
            mismatchKmer += 1
        else:
            # If the positions match then go to next kmer and next target
            # localisation
            if locaParsed[indParsed] == targetLoca[indTarg]:
                # print("Found !")
                indKmer += 1
                indTarg += 1
                indParsed = 0
            else:
                indParsed += 1
    return mismatchKmer


def link_kmer_fast(kmerList, locaList):
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
    kmerNb = len(kmerList)
    kmerLen = len(kmerList[0])
    fstLoca = locaList[0]
    lstLoca = locaList[-1]
    fstInd = 0
    lstInd = 0
    read = "".join(kmerList)
    # If first or last kmer have 0 localisation, return 0 localisation
    # for the read
    if len(lstLoca) == 1 and lstLoca[0] < 1:
        # print("No kmer localisation")
        return read, np.empty(1)
    if len(fstLoca) == 1 and fstLoca[0]:
        # print("No kmer localisation")
        return read, np.empty(1)
    locaRead = []
    while fstInd != len(fstLoca):
        if lstInd == len(lstLoca):
            # print("Next first localisation")
            fstInd += 1
            lstInd = 0
        else:
            expLoca = fstLoca[fstInd] + kmerLen * (kmerNb - 1)
            if expLoca >= lstLoca[lstInd]:
                if expLoca == lstLoca[lstInd]:
                    # print("Equal")
                    locaRead.append(fstLoca[fstInd])
                    fstInd += 1
                    lstInd = 0
                else:
                    # print("Superior")
                    lstInd += 1
            else:
                # wprint("Inferior")
                fstInd += 1
                lstInd = 0
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
