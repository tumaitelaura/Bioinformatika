from Bio import SeqIO
from Bio.Seq import Seq
from numpy import around
start_codons = ['ATG']
stop_codons = ['TAG', 'TAA', 'TGA']


def readFasta(handle: str):
    for seq_record in SeqIO.parse(handle, "fasta"):
        return seq_record


def getFrame(seq: Seq, i: int):
    return seq[i::3]


def lookForStartToStop(frame, stride, x, length):
    found = False
    for j in range(len(frame[stride:])):
        x += frame[j + stride]
        length += 1
        if frame[j + stride:j + stride + 3] in stop_codons:
            x += frame[j + stride + 1] + frame[j + stride + 2]
            found = True
            length += 2
            break
    return found, x


def lookForStopToStart(frame, stride, i, length):
    found = False
    x = ""
    for j in range(len(frame[stride:])):
        if frame[j + stride: j + stride + 3] in start_codons:
            x = frame[i: j + stride + 3]
            length += 3
            found = True
            j += 3
        elif frame[j + stride: j + stride + 3] in stop_codons:
            break
        length += 1
    return found, x


def findCodons(arr: list, frame: Seq, inArr: list, startOrStop: bool):
    for i in range(len(frame)):
        stride = i + 3
        length = 0
        if frame[i:stride] in inArr:
            length += 3
            if startOrStop:
                found, x = lookForStartToStop(
                    frame, stride, frame[i:stride], length)
            else:
                found, x = lookForStopToStart(frame, stride, i, length)
            if found:
                arr.append(x)
        i += length


def callThreeFrames(arr: list, seq: Seq, startOrStop):
    for i in range(3):
        frame = getFrame(seq, i)
        if startOrStop:
            findCodons(arr, frame, start_codons, startOrStop)
        else:
            findCodons(arr, frame, stop_codons, startOrStop)


def findCodonPairs(seq: Seq, startOrStop):
    arr = []
    callThreeFrames(arr, seq, startOrStop)
    rcSeq = seq.reverse_complement()
    callThreeFrames(arr, rcSeq, startOrStop)
    return arr


def filterByLength(arr: list):
    return [k for k in arr if len(k) >= 100]


def getAllPossible(seq: Seq, n1: int, n2: int):
    arr = []
    i = 0
    while i < len(seq):
        if seq[i:i+n2] not in arr:
            if len(seq[i:i+n2]) == n2:
                arr.append(seq[i:i+n2])
        i += n1
    return arr


def concatSeqs(seqs):
    sum = Seq("")
    for item in seqs:
        sum += item
    return sum


def getFreqs(seq: Seq, n1: int, n2: int):
    freqs = []
    allPossible = getAllPossible(seq, n1, n2)
    for possible in allPossible:
        freqs.append(Frequency(possible, seq.count(possible) / len(seq)))
    return freqs


def getGenomes(genomes: int):
    arr = []
    type = "bacterial"
    for i in range(genomes):
        num = i + 1
        if num > genomes / 2:
            type = "mamalian"
            num -= int(genomes / 2)
        arr.append(readFasta("data/" + type + str(num) + ".fasta"))
    return arr


def getPairs(genomes, startOrStop):
    arr = []
    for i in range(len(genomes)):
        arr.append(findCodonPairs(genomes[i].seq, startOrStop))
    return arr


class Frequency:
    def __init__(self, code, freq) -> None:
        self.code = code
        self.freq = freq


def getFreqsArr(genomes, codes1, codes2, n1, n2):
    arr = []
    for i in range(len(genomes)):
        concated = concatSeqs(codes1[i]) + concatSeqs(codes2[i])
        if concated:
            arr.append(getFreqs(concated, n1, n2))
        else:
            arr.append([])
    return arr


def findSameFreq(code, arr: list):
    for x in arr:
        if x.code == code:
            return x


def compareFreqs(arr1: list, arr2: list):
    resultArr = []
    for i in range(len(arr1)):
        sameCodeItem = findSameFreq(arr1[i].code, arr2)
        if sameCodeItem:
            delta = arr1[i].freq - sameCodeItem.freq
            if delta < 0:
                delta = delta * -1
            resultArr.append(delta)
    if len(resultArr) > 0:
        return (sum(resultArr) / len(resultArr)) * 1000
    return -1


def createPhilypMatrix(genomes, freqs):
    philypMatrix = []
    for i in range(len(genomes)):
        philypMatrix.append([])
    for i in range(len(genomes)):
        for j in range(len(genomes)):
            philypMatrix[i].append(compareFreqs(
                freqs[i], freqs[j]))
    return philypMatrix


def printPhilypMatrix(genomes, matrix):
    print(len(genomes))
    for i in range(len(genomes)):
        newList = around(matrix[i], 4)
        print(genomes[i].id, *newList)


# setup

genomes = getGenomes(8)

# 1.

startToStop = getPairs(genomes, True)

# 2.

stopToStart = getPairs(genomes, False)

# 3.

for i in range(len(genomes)):
    startToStop[i] = filterByLength(startToStop[i])
    stopToStart[i] = filterByLength(stopToStart[i])

# 4.

codonFreqsArr = getFreqsArr(genomes, startToStop, stopToStart, 1, 3)
dicodonFreqsArr = getFreqsArr(genomes, startToStop, stopToStart, 3, 6)

# 5.

printPhilypMatrix(genomes, createPhilypMatrix(genomes, codonFreqsArr))
printPhilypMatrix(genomes, createPhilypMatrix(genomes, dicodonFreqsArr))