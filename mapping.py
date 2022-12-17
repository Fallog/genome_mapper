import numpy as np
from Bio import SeqIO
from tqdm import tqdm
from Chromosome import Chromosome


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
    for base in read[
        ::-1
    ]:  # Do ... while, because have to check after the iteration if its good
        i_first_occ = chromo.first_occ[base]
        i_rank_mat = chromo.rank_mat[base]
        read = read[:-1]
        top = i_first_occ + i_rank_mat[top]
        try:
            if top < bottom:
                bottom = (
                    i_first_occ + i_rank_mat[bottom + 1] - 1
                )  # Because it will count 2 times the first
            else:
                bottom = i_first_occ + i_rank_mat[bottom]
        except IndexError:  # That also a way to eliminate wrong searching
            return np.array([-1])
        if bottom < top:
            return np.array([-1])
    return np.sort(chromo.suffix_table[top : bottom + 1], axis=0, kind="mergesort")


def link_kread(loc_list: np.ndarray, pat_len, nb_error=50000) -> np.ndarray:
    """Algorithm to link all kmer (kread here) depending on the number of error allowed

    Args:
        loc_list (np.ndarray): list of all localization of kmer
        pat_len (int): length on the kmer
        nb_error (int, optional): Number of error allowed. Defaults to 50000.

    Returns:
        np.ndarray: all the true localization depending of the number of error allowed
    """
    gap_start_end = 100 - pat_len
    km_restant = (100 // pat_len) - 2
    nb_bad_kmer = nb_error
    nb_good_needed = km_restant - nb_bad_kmer
    nb_good_kmer = 0
    res = []
    kpos = 1
    for i in loc_list[0]:
        for j in loc_list[-kpos]:
            if (j - i) == gap_start_end:
                if nb_bad_kmer >= km_restant:  # = max kmer ici
                    res.append(i)
                else:
                    nb_bad_kmer = nb_error
                    nb_good_needed = km_restant - nb_bad_kmer
                    nb_good_kmer = 0
                    res = []
                    kpos = 1
                    while nb_bad_kmer >= 0:
                        kpos += 1
                        gap_start_end -= pat_len
                        if check_link(i, loc_list, kpos, gap_start_end):
                            nb_good_kmer += 1
                        else:
                            nb_bad_kmer -= 1
                        if nb_good_kmer == nb_good_needed:
                            res.append(i)
                            break
    return res


def check_link(start, loc_list, k_pos, gap):
    for i in loc_list[-k_pos]:
        if (i - start) == gap:
            return True
    return False


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
