import numpy as np


def str_to_base(dna_seq: str) -> np.ndarray:
    """Transforms a dna_seq into an array of integers.

    Args:
        text (str): big string made only with A, T, C and G character

    Returns:
        np.ndarray: array resulting of the transformation
    """
    dna_seq = dna_seq.upper()
    res_array = np.empty(len(dna_seq), dtype=int)
    for i in range(len(dna_seq)):
        if dna_seq[i] == "$":
            res_array[i] = 1
        if dna_seq[i] == "A":
            res_array[i] = 2
        elif dna_seq[i] == "C":
            res_array[i] = 3
        elif dna_seq[i] == "G":
            res_array[i] = 4
        elif dna_seq[i] == "T":
            res_array[i] = 5

    return res_array


def inverse_sequence(dna_seq: str) -> str:
    """Returns the inversed complementary strand of the dnaSeq argument.
    'A' becomes 'T', 'C' becomes 'G' and inversely.

    Args:
        dna_seq (str): str made of the characters 'A', 'T', 'C' and 'G'

    Returns:
        str: complementary strand of dnaSeq, built backwards
    """
    seqLen = len(dna_seq)  # Performance
    dnaSeqInv = [0] * seqLen  # Performance

    for i in range(seqLen):
        base = dna_seq[i].capitalize()
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

    strand = "GCTTAGGAACTATACAGTT"
    invStrand = "AACTGTATAGTTCCTAAGC"
    print(f"Inversing test: {inverse_sequence(strand) == invStrand}")  # OK
