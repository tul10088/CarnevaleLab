import timeit
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import OPTICS

startTime = timeit.default_timer()
matrixHammings = np.load("/Users/aurelioaquila/PycharmProjects/Practice/hammings1000.npy")
print("1")
# Define the distance metric
def distance_metric(x, y):
    return np.sqrt(np.sum((x-y)**2))
print("2")

# Create the reachability graph
clustering = OPTICS(min_samples=10, max_eps=0.5, metric=distance_metric).fit(matrixHammings)
print("3")

# Compute the ordering
ordering = clustering.ordering_
print("4")

# Extract the clusters
labels = clustering.labels_
print("5")

# Print the number of clusters and noise points
num_clusters = len(set(labels)) - (1 if -1 in labels else 0)
num_noise_points = list(labels).count(-1)
print(f"Number of clusters: {num_clusters}")
print(f"Number of noise points: {num_noise_points}")

# Create a scatter plot of the data, colored by cluster
plt.scatter(matrixHammings[:, 0], matrixHammings[:, 1], c=labels, cmap='viridis')

# Add a colorbar to the plot
plt.colorbar()

# Add axis labels and a title to the plot
plt.xlabel('Feature 1')
plt.ylabel('Feature 2')
plt.title('OPTICS Clustering')

print("This algorithm takes: " + str(timeit.default_timer() - startTime) + " seconds")
