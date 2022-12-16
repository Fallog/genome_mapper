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
    if base == "$":
        return 0
    elif base == "A":
        return 1
    elif base == "C":
        return rank_mat["A"][-1] + 1
    elif base == "G":
        return rank_mat["A"][-1] + rank_mat["C"][-1] + 1
    else:
        return rank_mat["A"][-1] + rank_mat["C"][-1] + rank_mat["G"][-1] + 1


def string_search(read: str, chromo: Chromosome) -> np.ndarray:
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
    bottom = chromo.length - 1
    f_read = read
    for base in read[
        ::-1
    ]:  # Do ... while, because have to check after the iteration if its good
        i_first_occ = chromo.first_occ[base]
        i_rank_mat = chromo.rank_mat[base]
        read = read[:-1]
        top = i_first_occ + i_rank_mat[top]
        if top < bottom:
            bottom = (
                i_first_occ + i_rank_mat[bottom + 1] - 1
            )  # Because it will count 2 times the first
        else:
            bottom = i_first_occ + i_rank_mat[bottom]
        if not (top <= bottom):
            return np.array([-1])
    return np.sort(chromo.suffix_table[top : bottom + 1], axis=0, kind="mergesort")


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
    return [read[i : i + patt_len] for i in range(0, len(read), patt_len)]


def get_read_quality(read_loc: np.int64, loc_list: np.ndarray, patt_len: int) -> int:
    """Returns the "quality" of a read, i.e the number of localisation
    it is possible to retrieve between the 2 extreme localisations
    (first and last indexes of loc_list)


    Args:
        read_loc (np.int64): localisation of the read to be checked
        loc_list (np.ndarray): array of every localisation for each
            pattern of a given read
        patt_len (int): length of the pattern a read is divided
            into

    Returns:
        int: number of mismatched patterns in between first and last
            patterns in the patterns list
    """
    # Localisation where the kmer are supposed to be found
    targ_loc = [read_loc + (i * patt_len) for i in range(1, patt_len)]

    # Counter giving the number of unplaced kmer between the extreme reads
    nb_mismatch = 0
    index_patt = 1  # Index specifying on which kmer we are
    index_parsed = 0  # Index to parse the localisations of inside loc_list
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


def link_kmer_fast(patt_list: list[str], loc_lis: np.ndarray) -> np.ndarray:
    """Returns the reconstructed read and its localisation on a
    chromosome.
    Linkage of the 2 extreme kmers of kmerList (first and last) consists
    in finding if the localisation of the last is equal to the
    localisation of the first + N-1 times the length of a kmer, (N being
    the number of kmer in kmerList).

    Args:
        kmerlist (list[str]): list of every kmer of a given read, output
            of cut_read_to_kmer function
        locaList (np.ndarray): sorted array of every localisation for each
            kmer in kmerList, output of string_search function

    Returns:
        ndarray: localisation of the read on a given chromosome,
            np.empty(1) if no localisation is found
    """
    patt_nb = len(patt_list)
    patt_len = len(patt_list[0])
    # Every localisations of the first pattern
    first_loc = loc_lis[0]
    # Every localisation of the first pattern
    last_loc = loc_lis[-1]
    # Index to parse localisations of the first pattern
    first_index = 0
    # Index to parse localisations of the first pattern
    last_index = 0
    # Read is the concatenation of all the elements inside pattern list
    read = "".join(patt_list)
    # If first or last kmer have 0 localisation, return 0 localisation
    # for the read
    if len(last_loc) == 1 and last_loc[0] < 1:
        # print("No kmer localisation")
        return read, np.array([-1])
    if len(first_loc) == 1 and first_loc[0]:
        # print("No kmer localisation")
        return read, np.array([-1])
    locaRead = []
    while first_index != len(first_loc):
        if last_index == len(last_loc):
            # print("Next first localisation")
            first_index += 1
            last_index = 0
        else:
            exp_loc = first_loc[first_index] + patt_len * (patt_nb - 1)
            if exp_loc >= last_loc[last_index]:
                if exp_loc == last_loc[last_index]:
                    # print("Equal")
                    locaRead.append(first_loc[first_index])
                    first_index += 1
                    last_index = 0
                else:
                    # print("Superior")
                    last_index += 1
            else:
                # print("Inferior")
                first_index += 1
                last_index = 0
    return np.array(locaRead)


if __name__ == "__main__":
    chromo = []
    with tqdm(total=15, desc="Chromosomes importation") as pbar:
        for record in SeqIO.parse("SEQUENCES/P_fal_genome.fna", format="fasta"):
            chromo.append(str(record.seq))
            pbar.update(1)
    chromo1 = Chromosome("P_fal_chromosome_1", chromo[0], 1)

    # First read of mapping_P_fal.sam file
    test_read = "AAACCCTGAACCCTAAACCCTGAACCCTAAACCCTAAACCCTGAACCCTAAACCCTAAACCCTGAACCCTAAACCCTGAAACCTAAAACCTGAACCCTAA"
    kmer_fst_read = cut_read_to_kmer(test_read, 10)
    print(f"Kmers: {kmer_fst_read}")

    # testing for the first chromosome
    bwt_chr1 = chromo1.bwt
    # print(f"bwt: {bwt_chr1[:10000]}, type: {type(bwt_chr1)}")

    rank_mat = chromo1.rank_mat
    # print(f"First occ of C {get_first_occ(rank_mat, 'C')}")  # OK
    chromo1.compute_first_occurency()
    print(chromo1.first_occ)
    suffix_t = chromo1.suffix_table

    # Localisations of the first kmer
    locs = string_search(kmer_fst_read[0], chromo1)
    print(f"Locs first kmer: {locs}")
    verification_pattern(chromo1.DNA_dol, kmer_fst_read[0], locs)

    # Localisation of all the kmers
    loc_kmer = []
    for kmer in kmer_fst_read:
        locs = string_search(kmer, chromo1)
        loc_kmer.append(locs)
    print(f"Locs read: {loc_kmer}")

    # Read reconstruction
    loc_read = link_kmer_fast(kmer_fst_read, loc_kmer)
    print(f"Matching localisation(s) for the read: {loc_read}")

    print(type(loc_read[0]))
    # Read quality
    read_qlty1 = get_read_quality(loc_read[0], loc_kmer, patt_len=10)
    read_qlty2 = get_read_quality(loc_read[1], loc_kmer, patt_len=10)
    print(f"Quality of first localisation: {read_qlty1}")  # 1
    print(f"Quality of second localisation: {read_qlty2}")  # 5
