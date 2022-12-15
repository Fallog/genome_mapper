import numpy as np
from Bio import SeqIO
# import cProfile
from tqdm import tqdm
# import bwt
from Chromosome import Chromosome


def verification_pattern(chromo, kmer, locs):
    kmerL = kmer.lower()
    print(
        f"Number of kmer in the chromosome: {chromo.count(kmerL)}")
    print(f"Kmer length: {len(kmerL)}")
    i = 1
    for loc in locs:

        print(f"{i:>3} loc: {loc}")
        print(
            f"{chromo[loc -10:loc]}--{chromo[loc:loc + len(kmer)]}--{chromo[loc + len(kmer) +1:loc + len(kmer) + 11]}")
        i += 1


def compare_ndarray(array1, array2):
    resArray = []
    length = len(array1)
    for i in range(length):
        resArray.append(array1[i] == array2[i])
    return resArray


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
        ndarray: every localisation where kmer can be read. If no
        localisation is found, it is an "empty" ndarray -> np.empty(1)
    """
    bwtSort = sorted(list(bwtDna))
    lenBwt = len(bwtDna)
    # print(f"Sorted BWT: {bwtSort}")

    firstBase = bwtDna[0]
    # print(f"First letter: {firstBase}")

    bottom = 0
    top = lenBwt - 1  # Stay in the string boundaries
    i = len(kmer) - 1  # Stay in the kmer boundaries
    # print(f"first bottom: {bottom},  first top: {top}")

    nbOccur = 1

    while nbOccur > 0 and i >= 0:
        X = kmer[i]

        if X not in bwtDna:
            return kmer, np.empty(1)
        else:
            # print(f"i: {i}, X: {X}")

            firstOcc = bwtSort.index(X) - 1  # $ not counted
            # print(f"First index of {X} in bwtSort: {firstOcc}")

            # Counts all the occurences of X in bwtDna between an empty
            # rank and the following rank to deduce the rank of an
            # actual index e
            preBot = bottom
            rankBot = rankMat[X][preBot]
            # print(f"rankBot before while: {rankBot}")
            countBot = 0  # number of times the letter X is met in bwtDna
            # We scan bwtDna for X letter until we reach
            while rankBot == -1:
                preBot += 1  # e < f, we increase to stay in the boundaries
                if bwtDna[preBot] == X:
                    countBot += 1
                rankBot = rankMat[X][preBot]
            rankBot -= countBot
            # print(f"rankBot after while: {rankBot}")

            preTop = top
            rankTop = rankMat[X][preTop]
            # print(f"rankTop before while: {rankTop}")
            countTop = 0
            while rankTop == -1:
                preTop -= 1  # f > e, we decrease to stay in the boundaries
                if bwtDna[preTop] == X:
                    countTop += 1
                rankTop = rankMat[X][preTop]
            rankTop += countTop
            # print(f"rankTop after while: {rankTop}")

            # The first character of the BWT has a rank of 1 in the
            # rank matrix but it is the first appearance of this
            # character. To take into account this information, we
            # decrease e of 1 in presence of this character, i.e when
            # there is no character before it
            if bottom == 0 and X == firstBase:
                bottom = firstOcc + rankBot - 1
            else:
                bottom = firstOcc + rankBot
            top = firstOcc + rankTop
            # print(f"bottom: {bottom}, top: {top}")

            nbOccur = top - bottom  # quantity of elements between 2 indexes

            if nbOccur == 0:
                # print("Kmer not found !")
                return kmer, np.empty(1)
            # print(f"Number of pattern: {nbOccur}\n")

            i -= 1

            # True if the entire pattern is crossed
            if i == -1:
                # Positions of the kmer in the chromosome computed only when
                # the kmer totally localised to save time
                locs = suffixTab[bottom:top]
                return kmer, np.sort(locs, kind="mergesort")

    return kmer


def cut_read_to_kmer(read: str, kLen: int):
    """Divide the given read argument into a list of k-mer,
    smaller strings of size the k_len argument.

    Args:
        read (str): str made of the characters 'A', 'T', 'C' and 'G'
        kLen (int): size of the k-mer to be created from the read
            sequence

    Returns:
        list[str]: list of all the k-mer created from the read
    """
    readLen = len(read)  # performance
    kmerList = [0] * (readLen // kLen)
    readCnt = 0
    kCnt = 0
    while readCnt <= readLen - kLen:
        kmerList[kCnt] = (read[readCnt:readCnt + kLen])
        readCnt += kLen
        kCnt += 1

    # Adding remaining nucleotides in case of non divisible k_len
    if readLen % kLen == 0:
        return kmerList
    else:
        return kmerList + read[readCnt:]


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


def get_read_quality(locaFst, locaList, kmerLen) -> int:
    """_summary_

    Args:
        locaFst (_type_): _description_
        locaList (_type_): _description_
        kmerLen (_type_): _description_

    Returns:
        int: _description_
    """
    targetLoca = [locaFst + (i * kmerLen) for i in range(1, kmerLen)]
    print(f"Locs we are looking for: {targetLoca}")
    # Counter giving the number of unplaced kmer between the extreme reads
    mismatchKmer = 0
    indKmer = 1  # Index specifying on which kmer we are
    indParsed = 0  # Index to parse the localisations of inside locaList
    indTarg = 0  # Index to parse targetLoca
    while indKmer < kmerLen - 1:
        locaParsed = locaList[indKmer]
        print(f"Locations parsed: {locaParsed} len: {len(locaParsed)}")
        # If it reaches the end of actual localisation list
        # before returning anything go to next kmer but count one
        # unplaced kmer
        if indParsed == len(locaParsed):
            print(f"indParsed: {indParsed}")
            indKmer += 1
            indTarg += 1
            indParsed = 0
            mismatchKmer += 1
        else:
            # If
            if locaParsed[indParsed] == targetLoca[indTarg]:
                print("Found !")
                indKmer += 1
                indTarg += 1
                indParsed = 0
            else:
                indParsed += 1
    return mismatchKmer


def link_kmer_fast(kmerlist: list[str], locaList, maxError=0):
    kmerNb = len(kmerlist)
    kmerLen = len(kmerlist[0])
    fstLoca = locaList[0]
    lstLoca = locaList[-1]
    fstInd = 0
    lstInd = 0
    read = "".join(kmerlist)
    # If first or last kmer have 0 localisation, return 0 localisation
    # for the read
    if len(lstLoca) == 1 and lstLoca[0] < 1:
        print("No kmer localisation")
        return read, np.empty(1)
    if len(fstLoca) == 1 and fstLoca[0]:
        print("No kmer localisation")
        return read, np.empty(1)
    locaRead = []
    while fstInd != len(fstLoca):
        if lstInd == len(lstLoca):
            print("Next first localisation")
            fstInd += 1
            lstInd = 0
        else:
            expLocaMin = fstLoca[fstInd] + kmerLen * (kmerNb - 1) - maxError
            expLocaMax = expLocaMin + (maxError * 2)
            if expLocaMin >= lstLoca[lstInd] or expLocaMax >= lstLoca[lstInd]:
                if expLocaMin == lstLoca[lstInd] or expLocaMax == lstLoca[lstInd]:
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

    reads = []
    with tqdm(total=1500000, desc="Reads importation") as pbar:
        for record in tqdm(SeqIO.parse("SEQUENCES/P_fal_reads.fq", format="fastq")):
            reads.append((record.seq, record.id))
            pbar.update(1)

    # First read of mapping_P_fal.sam file
    readTest = str(reads[100000 - 21515][0])
    print(f"First read: {readTest} type: {type(readTest)}")
    kmerFstRead = cut_read_to_kmer(readTest, 10)
    print(f"First kmer: {kmerFstRead}")

    # testing for the first chromosome
    bwtChromo1 = chromo1.bwt
    # print(f"bwt: {bwtChromo1[:10000]}, type: {type(bwtChromo1)}")

    bwtList = list(bwtChromo1)
    # print(f"Sorted BWT list: {sorted(bwtList)[:100]}")

    # cProfile.run("rankMat = bwt.create_rank_mat(bwtChromo1)")
    rankMat = chromo1.rank_mat.item()
    # print(f"Rank matrix of A: {rankMat['A'][260000:261000]}")

    locs = []
    for kmer in kmerFstRead:
        locs.append(search_kmer_pos(bwtChromo1, rankMat,
                    chromo1.suffix_table, kmer)[1])
    print(f"Kmer locs on chromo1: {locs} length: {len(locs)}")

    # verification_pattern(chromo1.DNA, kmerFstRead[0], locs[0])
    # recoRead = link_kmer(kmerFstRead, locs)
    recoReadFst = link_kmer_fast(kmerFstRead, locs)
    print(f"Reconstructed read: {recoReadFst}")
    qltyPos1 = get_read_quality(recoReadFst[1][0], locs, 10)  # 1
    qltyPos2 = get_read_quality(recoReadFst[1][1], locs, 10)  # 5
    print(f"Read qlty1: {qltyPos1} 2: {qltyPos2}")
    # print(f"Reconstructed read: {recoRead} Lenght read: {len(recoRead[0])}")
    # print(f"Actual read: {readTest:>7}")
