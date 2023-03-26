import timeit
import Bio.AlignIO  #reading fasta file
import numpy as np  #data processing
import pandas as pd #data processing
import matplotlib.pyplot as plt #visualization of data
from sklearn.preprocessing import OneHotEncoder #one hot encoding of amino acids
from sklearn.cluster import DBSCAN  #clustering method
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
        #print(record.seq)
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
    np.save("./NEWhammings100.npy", matrixHammings)

    # print(matrixHammings.flatten())
    # print(len(matrixHammings))

    # counts, bins = np.histogram(matrixHammings, bins= 100)
    #
    #
    # # print(matrixHammings)
    #
    # plt.stairs(counts, bins)
    # plt.show()


    # Cluster the dataset with OPTICS
    #def distance_metric(x, y):
        #return np.sqrt(np.sum((x-y)**2))
    #clustering = OPTICS(min_samples=5, max_eps=120, metric=distance_metric).fit(matrixHammings)
    clustering = DBSCAN(eps=120, min_samples=5, metric='precomputed').fit(matrixHammings)
    labels = clustering.labels_
    print(clustering)
    print(clustering.labels_)


    # Visualize the results
    #plt.scatter(matrixHammings[:, 0], matrixHammings[:, 1], c=labels, cmap='viridis')
    #plt.show()


# if hasattr(clustering, "labels_"):
    #     y_pred = clustering.labels_.astype(int)
    # else:
    #     y_pred = clustering.predict(matrixHammings)
    #
    # print(y_pred)

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    print("Estimated number of clusters: %d" % n_clusters_)
    print("Estimated number of noise points: %d" % n_noise_)




    #Uniprot to GO (tells you voltage gated postassium channel); Find FTP server with uniprot to GO via google (GO = Gene Ontology)
    #One cluster for every Uniprot code (the name of each sequence in the FASTA file ex. A0A2R2MPD5 for seq 1, 2, 3... etc.)

    # count = 0
    #
    # for i in clustering.labels_:
    #     if (i == -1):
    #         count+=1
    #
    # print(len(np.unique(clustering.labels_)))
    # print(count)

    # newF = open("testFileWrite.txt", 'w')
    # newF.write(str(matrixHammings))

def getIDs():
    newF = open(r'/Users/aurelioaquila/Documents/America/TEMPLE/Research/PF00520_rp15.txt')
    fileIDset = set()
    fileIDs = open("sequenceIDs.txt",'w')
    for line in newF:
        if '>' in line:
            seqRE = re.findall(r"[^>/]+", line)
            fileIDset.add((seqRE[0]) + "\n")
    set(seqRE)
    #print(fileIDset)
    for line in fileIDset:
        fileIDs.write(line)

startTime = timeit.default_timer()
main()
print("This algorithm takes: " + str(timeit.default_timer() - startTime) + " seconds")


def test():
    newF = open("testFileWrite.txt", 'w')
    a = np.array([[2, 3], [4, 5]])
    b = np.transpose(a)
    c = np.dot(a, b)
    newF.write(str(c))

#think about options to save the data in a txt file or something else
#next step is analysis of data and clustering, after optimization of the algorithm