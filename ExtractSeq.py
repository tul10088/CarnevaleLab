import timeit
import pandas as pd

# read the Excel file into a pandas dataframe
df = pd.read_excel('/Users/aurelioaquila/Documents/America/TEMPLE/Research/uniprot.xlsx')

# create a list of tuples from the dataframe, where each tuple represents a row of data
data_list = [(row['Entry'], row['Protein names']) for index, row in df.iterrows()]

# print the list
print(data_list)

groups = []
uniquekeys = []
for k, g in groupby(data, lambda x: x[1]):
    groups.append(list(g))      # Store group iterator as a list
    uniquekeys.append(k)

#print(groups)

startTime = timeit.default_timer()
print("This algorithm takes: " + str(timeit.default_timer() - startTime) + " seconds")
