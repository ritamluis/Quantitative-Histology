import os

import numpy as np

import skimage as ski
from skimage import io, filter, color, exposure, morphology, feature, measure, transform

from sklearn.linear_model import BayesianRidge
from sklearn.mixture import GMM

from scipy import spatial, ndimage, signal, stats

import matplotlib.pyplot as plt




dirq = "data/"
files = []
for i in os.listdir(dirq) :
    if i.endswith(".jpg") :
        files.append(dirq + i)


for f in files :

	image = exposure.equalize_adapthist(io.imread(f))

	binary = filter.threshold_adaptive(exposure.adjust_sigmoid(image[:, :, 0], cutoff=0.4, gain = 30), 301).astype(bool)
	clean = morphology.binary_closing(binary, morphology.disk(3)).astype(bool)
	clean = morphology.remove_small_objects(clean, 200)
	clean = morphology.remove_small_objects( (1-clean).astype(bool), 200)

	local_density = filter.gaussian_filter(clean, 61)

	ent = filter.gaussian_filter(filter.rank.entropy(local_density, morphology.disk(3)), 75)

	ent -= ent.min()
	ent /= ent.max()

	local_density -= local_density.min()
	local_density /= local_density.max()

	info = ent * (1 + local_density)

	bw = (info) > filter.threshold_otsu(info)


	C = measure.find_contours(bw, 0.5)
	centroid = []
	vals = []

	for c in C :
	    centroid.append(np.linalg.norm([c[:, 1].mean() - bw.shape[1] / 2, c[:, 0].mean() - bw.shape[0] / 2]))
	    vals.append(local_density.T[c.astype(int)].sum())

	cent = C[np.argmin(centroid / np.array(vals))]


	fix, ax = plt.subplots()
	plt.imshow(image)
	ax.plot(cent[:, 1], cent[:, 0], lw=5, c="k", alpha = 0.7)

	plt.savefig(f.split("Sheep")[1])
	print ("Done " + f.split("Sheep")[1])
