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
        list[str, int]: list of all the k-mer created from the read
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


def link_kmer(kmerList):
    read = ""
    for kmerTuple in kmerList:
        kmerOrder = kmerTuple[2]
    return read


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
    fst_read = """TTTCCTTTTTAAGCGTTTTATTTTTTAATAAAAAAAATATAGTATTATATAGTAACGGGTGAAAAGATCCATATAAATAAATATATGAGGAATATATTAA"""
    print(f"Cuttinng test: {cut_read_to_kmer(fst_read, 20, 1)}")  # OK

    frag = "TTTCCTTTTT"
    invFrag = "AAAAAGGAAA"
    print(f"Inversing test: {inverse_sequence(frag) == invFrag}")  # OK
    strand = "GCTTAGGAACTATACAGTT"
    invStrand = "AACTGTATAGTTCCTAAGC"
    print(f"Inversing test: {inverse_sequence(strand) == invStrand}")  # OK
