from collections import OrderedDict
import pandas as pd
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster
from sklearn.cluster import AgglomerativeClustering
from sklearn import manifold
from time import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter
import umap.umap_ as umap

# set pandas options for ease of use
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

# read in the data (note that we are reading in a file which already has the exact columns we want)
df = pd.read_csv('bulbar_5.csv')

# drop rows with missing values (not that I believe there are any, just in case)
df = df.dropna()

# reset the index
df = df.reset_index(drop=True)
print(df)

# split off numerical data
df_num = df.loc[:, (df.columns != 'Participant_ID') & (df.columns != 'SiteOnset_Class')]

# normalize all the gene expression data using the Standard Scaler
scaler = StandardScaler()
scaler.fit(df_num)
scaled = scaler.fit_transform(df_num)
df_num = pd.DataFrame(scaled, columns=df_num.columns)
print(df_num)

# data is now ready for agglomerative clustering
clusterer = AgglomerativeClustering(distance_threshold=0, n_clusters=None)
clusters = clusterer.fit(df_num)

# clustering has now been performed and is ready for visualization

# this dendrogram-plotting function is adapted from a scikit learn tutorial
def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, orientation='left', **kwargs)
    return linkage_matrix


# plot the dendrogram
fig, ax = plt.subplots()
plt.title('Hierarchical Clustering Dendrogram')
Z = plot_dendrogram(clusters, truncate_mode='level', p=3, color_threshold=0.0)
plt.xlabel("Cophenetic Distance")
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
plt.ylabel("Index or size of datapoint in cluster")
plt.show()

k = 2
R = fcluster(Z, k, criterion='maxclust')

# look at the data in the clusters when there are exactly 2 clusters (since we hope they split into bulbar vs. limb)
df_num['cluster'] = R

# color code the data according to which cluster it is in
cluster_colors = [None] * len(df_num.index)
for index, row in df_num.iterrows():
    if row['cluster'] == 1:
        cluster_colors[index] = 'blue'
    elif row['cluster'] == 2:
        cluster_colors[index] = 'green'

# separately, color code the data based on whether it is bulbar or limb
# limb = 1, bulbar = 0
onset_colors = [None] * len(df.index)
for index, row in df.iterrows():
    if row['SiteOnset_Class'] == 0: #bulbar
        onset_colors[index] = 'blue'
    elif row['SiteOnset_Class'] == 1: #limb
        onset_colors[index] = 'green'

# perform t-SNE and plot results
Axes3D

n_points = 1000
X = df_num.loc[:, df_num.columns != 'cluster']
color = 'green'
n_neighbors = 10
n_components = 2

# Create figure
fig = plt.figure(figsize=(15, 8))
fig.suptitle("Dimensionality Reduction", fontsize=16)

# Set-up manifold methods for dimensionality reduction
methods = OrderedDict()
methods['MDS'] = manifold.MDS(n_components=n_components, max_iter=100, n_init=1)
methods['UMAP'] = umap.UMAP()
methods['t-SNE'] = manifold.TSNE(n_components=n_components, init='pca', random_state=0)

# Plot results
for i, (label, method) in enumerate(methods.items()):
    t0 = time()
    Y = method.fit_transform(X)
    t1 = time()
    print("%s: %.2g sec" % (label, t1 - t0))
    ax = fig.add_subplot(2, 5, 2 + i + (i > 3))
    ax.scatter(Y[:, 0], Y[:, 1], c=cluster_colors, cmap=plt.cm.Spectral)
    ax.set_title("%s (%.2g sec)" % (label, t1 - t0))
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.axis('tight')

plt.show()

# compare MDS clustering with actual bulbar vs. limb
Y = methods['MDS'].fit_transform(X)
fig, axs = plt.subplots(ncols=2)
axs[0].scatter(Y[:, 0], Y[:, 1], c=cluster_colors, cmap=plt.cm.Spectral)
axs[0].set_title("MDS")
axs[0].xaxis.set_major_formatter(NullFormatter())
axs[0].yaxis.set_major_formatter(NullFormatter())
axs[0].axis('tight')

axs[1].scatter(Y[:, 0], Y[:, 1], c=onset_colors, cmap=plt.cm.Spectral)
axs[1].set_title("Bulbar vs. Limb")
axs[1].xaxis.set_major_formatter(NullFormatter())
axs[1].yaxis.set_major_formatter(NullFormatter())
axs[1].axis('tight')

plt.show()

# compare UMAP clustering with actual bulbar vs. limb
Y = methods['UMAP'].fit_transform(X)
fig, axs = plt.subplots(ncols=2)
axs[0].scatter(Y[:, 0], Y[:, 1], c=cluster_colors, cmap=plt.cm.Spectral)
axs[0].set_title("UMAP")
axs[0].xaxis.set_major_formatter(NullFormatter())
axs[0].yaxis.set_major_formatter(NullFormatter())
axs[0].axis('tight')

axs[1].scatter(Y[:, 0], Y[:, 1], c=onset_colors, cmap=plt.cm.Spectral)
axs[1].set_title("Bulbar vs. Limb")
axs[1].xaxis.set_major_formatter(NullFormatter())
axs[1].yaxis.set_major_formatter(NullFormatter())
axs[1].axis('tight')

plt.show()

# compare t-SNE clustering with actual bulbar vs. limb
Y = methods['t-SNE'].fit_transform(X)
fig, axs = plt.subplots(ncols=2)
axs[0].scatter(Y[:, 0], Y[:, 1], c=cluster_colors, cmap=plt.cm.Spectral)
axs[0].set_title("t-SNE")
axs[0].xaxis.set_major_formatter(NullFormatter())
axs[0].yaxis.set_major_formatter(NullFormatter())
axs[0].axis('tight')

axs[1].scatter(Y[:, 0], Y[:, 1], c=onset_colors, cmap=plt.cm.Spectral)
axs[1].set_title("Bulbar vs. Limb")
axs[1].xaxis.set_major_formatter(NullFormatter())
axs[1].yaxis.set_major_formatter(NullFormatter())
axs[1].axis('tight')

plt.show()
