# import cProfile
# import random
import numpy as np


def strToAscii(text: str):  # USELESS
    """Transform a string to is Ascii value ant put it in a numpy array
    TO DO : update for DNA (A --> 1 : C --> 2 : G --> 3 : T --> 4)

    Args:
        text (str): DNA sequence

    Returns:
        array: Number array (representing the DNA sequence)
    """
    res_array = np.frombuffer(text.encode(), np.int8)
    return res_array


def strToBase(text: str):
    text = text.upper()
    res_array = np.empty(len(text), dtype=int)
    for i in range(len(text)):
        if text[i] == "$":
            res_array[i] = 1
        if text[i] == "A":
            res_array[i] = 2
        elif text[i] == "C":
            res_array[i] = 3
        elif text[i] == "G":
            res_array[i] = 4
        elif text[i] == "T":
            res_array[i] = 5

    return res_array


def asciiToStr(ascii_list) -> str:
    """Transform a Numpy array of int to a string
    TO DO : adapt to ACGT

    Args:
        ascii_list (array): Numpay array of int

    Returns:
        str: String representative of the initial array
    """
    res_txt = ascii_list.tobytes().decode("ascii")
    return res_txt


def printDc3Var(p12, r12, r12s, index, tp, iter=-1) -> None:
    param = locals()
    str_to_print = ""
    for key, val in param.items():
        if key == "iter":
            if iter >= 0:
                print(f"---- ITERATION N°{iter} ----\n")
        else:
            str_to_print += f"Value of {key} is actually : {val}\n"
    print(str_to_print)


def cut_read_to_kmer(read: str, kLen: int, readId):
    """Divide the given read argument into a list of k-mer, i.e. smaller strings,
    of size the k_len argument.

    Args:
        read (str): str made of the characters 'A', 'T', 'C' and 'G'
        kLen (int): size of the k-mer to be created from the read
            sequence
        readId (int): unique identifier specifying to which read the
            kmer belongs

    Returns:
        list[str, int, int]: list of all the k-mer created from the read
            along with readId and the order of the kmer on the read
    """
    readLen = len(read)  # performance
    kmerList = [0] * (readLen // kLen)
    readCnt = 0
    kCnt = 0
    while readCnt <= readLen - kLen:
        kmerList[kCnt] = (read[readCnt:readCnt + kLen], readId, kCnt)
        readCnt += kLen
        kCnt += 1

    # Adding remaining nucleotides in case of non divisible k_len
    return kmerList + [(read[readCnt:], readId, kCnt + 1)]


def link_kmer(kmerList, locaList):
    indFirst = 0  # index to parse the localisation of the first kmer
    indLoca = 0  # index to parse the localisation of the following kmers
    indKmer = 1  # index specifying on which following kmer we are
    locaFirst = locaList[0]  # localisations list of the first kmer
    locaNext = locaList[indKmer]  # localisations list of the following kmer
    read = kmerList[0]
    readLen = len(kmerList)  # number of kmer in a read
    kmerLen = len(kmerList[0])  # kmer have all the same size
    while indKmer != readLen:
        # Si on arrive au bout des localisations, on passe à la
        # localisation suivante du premier kmer en retournant à
        # la première localisation du kmer suivant (le 1 donc)
        if indLoca == len(locaNext):
            indFirst += 1
            indKmer = 1
            indLoca = 0
            read = kmerList[0]
        else:
            # Si les localisations se suivent, passer au kmer suivant
            # tout en revenant à sa première localisation
            if locaFirst[indFirst] + kmerLen * indKmer == locaNext[indLoca]:
                read += kmerList[indKmer]
                indKmer += 1
                indLoca = 0

            else:  # sinon regarder la localisation suivante du kmer actuel
                indLoca += 1
        locaNext = locaList[indKmer]
    return read, locaFirst[indFirst]


def inverse_sequence(dnaSeq: str) -> str:
    """Returns the inversed complementary strand of the dnaSeq argument.
    'A' becomes 'T', 'C' becomes 'G' and inversely.

    Args:
        dnaSeq (str): str made of the characters 'A', 'T', 'C' and 'G'

    Returns:
        str: complementary strand of dnaSeq, built backwards
    """
    seqLen = len(dnaSeq)  # Performance
    dnaSeqInv = [0] * seqLen  # Performance

    for i in range(seqLen):
        base = dnaSeq[i].capitalize()
        # dnaSeqInv is built backwards because it the inversed of dnaSeq
        if base == "A":
            dnaSeqInv[-(i + 1)] = "T"
        elif base == "T":
            dnaSeqInv[-(i + 1)] = "A"
        elif base == "C":
            dnaSeqInv[-(i + 1)] = "G"
        else:
            dnaSeqInv[-(i + 1)] = "C"
    return "".join(dnaSeqInv)  # Recreate a string


if __name__ == "__main__":
    testRead = """TTTCCTTTTTAAGCGTTTTATTTTTTAATAAAAAAAATATAGTATTATATAGTAACGGGTGAAAAGATCCATATAAATAAATATATGAGGAATATATTAA"""
    print(f"Cutting test: {cut_read_to_kmer(testRead, 20, 1)}")  # OK
    kmerList = cut_read_to_kmer(testRead, 10, 1)

    strand = "GCTTAGGAACTATACAGTT"
    invStrand = "AACTGTATAGTTCCTAAGC"
    print(f"Inversing test: {inverse_sequence(strand) == invStrand}")  # OK
