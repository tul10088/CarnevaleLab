import timeit
import Bio.AlignIO  #reading fasta file
import numpy as np  #data processing
import pandas as pd #data processing
import matplotlib.pyplot as plt #visualization of data
from sklearn.preprocessing import OneHotEncoder #one hot encoding of amino acids
from sklearn.cluster import OPTICS  #clustering method
from sklearn import metrics #counting number of clusters and outliers (noise points)
import seaborn as sns   #visualization of clustering

def main():
#This block of code opens the text file so that python can read and parse the protein sequences
    file = open(r'/Users/aurelioaquila/Documents/America/TEMPLE/Research/PF00520_rp15.txt') #original fasta file
    f = open("SeqPF00520_rp15.txt", 'w')
    fileWrite = open("SeqPF00520_rp15.txt", 'r') #file with all sequences line by line
    alignment = Bio.AlignIO.read(file, "fasta")
    listSeq = []

#This creates a list that has a sequence at each index, by this I mean that sequence 1 is the first position in the list, sequence 2 is the next position, and so on.
    for record in alignment:
        f.write(str(record.seq) + "\n")
        listSeq.append(str(record.seq).upper())
    #print(listSeq[0])

    aminoAcidList = ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W', 'Y', 'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q', 'X', '-']
#This block of code creates the binary arrays (vectors) for each character in the aminoAcidList
    df = pd.DataFrame(aminoAcidList)
    ohe = OneHotEncoder()
    arraySequences = ohe.fit_transform(df).toarray()
    #print(arraySequences)

#This dictionary attaches the amino acid character from aminoAcidList to the corresponding binary vector
    newDict = {}
    for i in range(len(aminoAcidList)):
        newDict.update({aminoAcidList[i] : arraySequences[i]})
    #print(newDict)

#This block of code currently appends each protein sequences' binary vector to a list called binaryList.
#Where index 0 is every char in seq1 and index 1 is every char in se2, etc.
#benchmark this code block = 3.26s
    binaryList = []
    tempList = []
    for item in listSeq:
        for character in item:
            tempList.append((newDict[character]))
        binaryList.append(tempList)
        tempList = []
    #print(binaryList[0])
    #print(len(binaryList))

#This block of code will compute the hamming distance via the matrices gotten from above by changing the 3D matrix to a 2D matrix via flattening
    matrixCharBySeq = np.asarray(binaryList)
    x = matrixCharBySeq
    print(x.shape)
    a= np.reshape(x, (20081, 1642*22))
    print(a.shape)

    matrixSeqByChar = np.transpose(a)
    matrixDot = np.dot(a, matrixSeqByChar)
    matrixHammings = len(listSeq[0]) - matrixDot
    np.save("./NEWhammings.npy", matrixHammings)

startTime = timeit.default_timer()
main()
print("This algorithm takes: " + str(timeit.default_timer() - startTime) + " seconds")