import timeit
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import OPTICS

startTime = timeit.default_timer()
matrixHammings = np.load("/Users/aurelioaquila/PycharmProjects/Practice/hammings1000.npy")

# Fit the OPTICS algorithm to the data
optics = OPTICS(min_samples=50, xi=0.05, min_cluster_size=0.1)
optics.fit(matrixHammings)

# Plot the clusters
reachability = optics.reachability_
labels = optics.labels_
colors = [plt.cm.Spectral(each)
          for each in np.linspace(0, 1, len(set(labels)))]

plt.figure(figsize=(10, 7))
for klass, color in zip(range(0, len(colors)), colors):
    Xk = matrixHammings[labels == klass]
    reachabilityk = reachability[labels == klass]
    plt.plot(Xk[:, 0], reachabilityk, 'o', markerfacecolor=tuple(color),
             markeredgecolor='k', markersize=6)

plt.title('OPTICS clustering')
plt.xlabel('Data points')
plt.ylabel('Reachability distance')
plt.show()

print("This algorithm takes: " + str(timeit.default_timer() - startTime) + " seconds")
