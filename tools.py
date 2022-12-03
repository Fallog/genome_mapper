import cProfile
import random
import numpy as np


def strToAscii(text: str):  # ya un numpy loadtxt qui existe
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
    res_array = np.empty(len(text))
    for i in range(len(text)):
        if text[i] == "A":
            res_array[i] = 1
        elif text[i] == "C":
            res_array[i] = 2
        elif text[i] == "G":
            res_array[i] = 3
        else:
            res_array[i] = 4
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


def cut_read_to_kmer(read: str, k_len) -> list[str]:
    """Divide the given read argument into a list of k-mer, i.e. smaller strings,
    of size the k_len argument.

    Args:
        read (str): str made of the characters 'A', 'T', 'C' and 'G'
        k_len (_type_): size of the k-mer to be created from the read sequence

    Returns:
        list[str]: list of all the k-mer created from the read argument
    """
    read_len = len(read)  # performance
    kmer_l = [0] * (read_len // k_len)
    read_cnt = 0
    k_cnt = 0
    while read_cnt <= read_len - k_len:
        kmer_l[k_cnt] = read[read_cnt:read_cnt + k_len]
        read_cnt += k_len
        k_cnt += 1

    # Adding remaining nucleotides in case of non divisible k_len
    return kmer_l + [read[read_cnt:]]


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
        base = dnaSeq[i]
        # dnaSeqInv is built backwards because it the inversed of dnaSeq
        if base == 'A':
            dnaSeqInv[-(i + 1)] = 'T'
        elif base == 'T':
            dnaSeqInv[-(i + 1)] = 'A'
        elif base == 'C':
            dnaSeqInv[-(i + 1)] = 'G'
        else:
            dnaSeqInv[-(i + 1)] = 'C'
    return "".join(dnaSeqInv)  # Recreate a string


if __name__ == "__main__":
    fst_read = """TTTCCTTTTTAAGCGTTTTATTTTTTAATAAAAAAAATATAGTATTATATAGTAACGGGTGAAAAGATCCATATAAATAAATATATGAGGAATATATTAA"""
    print(f"Cuttinng test: {cut_read_to_kmer(fst_read, 20)}")

    frag = "TTTCCTTTTT"
    invFrag = "AAAAAGGAAA"
    print(f"Inversing test: {inverse_sequence(frag) == invFrag}")  # OK

    # liste = np.array(
    #     [
    #         (5, 28, 1),
    #         (9, 111, 9),
    #         (9, 1, 5),
    #         (99, 1, 9),
    #         (2, 1, 89),
    #         (29, 19, 1),
    #     ],
    #     dtype=[("1", int), ("2", int), ("3", int)],
    # )

    # random.seed(12)
    # datatest = [
    #     (random.randint(1, 4), random.randint(1, 4), random.randint(1, 4))
    #     for i in range(500000)
    # ]

    # datatest = np.array(datatest, dtype=[("1", int), ("2", int), ("3", int)])

    # # print(datatest)

    # print(liste)
    # res = quickSort(liste)
    # print(liste)
    # print(res)
    # cProfile.run("quickSort(datatest)")
