# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os

import numpy as np

import skimage as ski
from skimage import io, filter, color, exposure, morphology, feature, draw, measure, transform

from scipy.spatial import distance
from scipy import ndimage

import matplotlib
%matplotlib inline
import matplotlib.pyplot as plt
import seaborn
matplotlib.rcParams["figure.figsize"] = (16, 10)

# <codecell>

# Load filenames
dirq = "/Users/qcaudron/repositories/Quantitative-Histology/10x/data/"
files = []
for i in os.listdir(dirq) :
    if i.endswith(".jpg") :
        files.append(dirq + i)
        

offset = 1000
lengthscale = 301
A = io.imread(files[4])


# Things to measure for each image
in_count = []
out_count = []
in_size = []
out_size = []
in_area_density = []
out_area_density = []
in_point_density = []
out_point_density = []
var_point_density = []
cov_point_density = []
var_area_density = []
cov_area_density = []
fs_area = []
pairwise = []


image = exposure.equalize_adapthist(A)
io.imshow(image)
plt.grid(False)

# <codecell>

# Process for nuclei

binary = filter.threshold_adaptive(exposure.adjust_sigmoid(A[:, :, 0], cutoff=0.4, gain = 30), 301).astype(bool)
clean = morphology.binary_closing(binary, morphology.disk(3)).astype(bool)
clean = morphology.remove_small_objects(clean, 200)
clean = morphology.remove_small_objects( (1-clean).astype(bool), 200)

io.imshow(clean)
plt.grid(False)

# <codecell>

# Find contour of inflammatory zone

local_density = filter.gaussian_filter(clean, 61)
local_density -= local_density.min()
local_density /= local_density.max()

ent = filter.gaussian_filter(filter.rank.entropy(local_density, morphology.disk(3)), 75)
ent -= ent.min()
ent /= ent.max()

info = ent * (1 + local_density)

bw = (info) > filter.threshold_otsu(info)

C = measure.find_contours(bw, 0.5)
centroid = []
vals = []

for c in C :
    centroid.append(np.linalg.norm([c[:, 1].mean() - bw.shape[1] / 2, c[:, 0].mean() - bw.shape[0] / 2]))
    vals.append(local_density.T[c.astype(int)].sum())

cent = C[np.argmin(centroid / np.array(vals))]
path = matplotlib.path.Path(cent)

io.imshow(image)
plt.plot(cent[:, 1], cent[:, 0], lw=5, c="k", alpha = 0.7)
plt.grid(False)

# <codecell>

# Segmentation

distance = ndimage.distance_transform_edt(clean)
peaks = feature.peak_local_max(distance, indices=False, labels = clean)
markers = ndimage.label(peaks)[0]
labels = morphology.watershed(-distance, markers, mask=clean)

io.imshow(labels, interpolation = "nearest", cmap = plt.cm.spectral)
plt.grid(False)

# <codecell>

# Properties of each nucleus
nuclei = measure.regionprops(labels)
contour = matplotlib.path.Path(cent, closed=True)

incount = 0
outcount = 0
insize = []
outsize = []
outer_nuclei = []

# Number of nuclei, and their sizes, inside and outside the focus
for n in nuclei :
    if contour.contains_point(n.centroid) :
        incount += 1
        insize.append((labels == n.label).sum())
    else :
        outcount += 1
        outsize.append((labels == n.label).sum())
        outer_nuclei.append(n)

in_count.append(incount)
out_count.append(outcount)
in_size.append(np.mean(insize))
out_size.append(np.mean(outsize))




# Using the shoelace algorithm, compute the area of the focus
x = path.vertices[:, 0]
y = path.vertices[:, 1]
x = np.append(x, path.vertices[0, 0])
y = np.append(y, path.vertices[0, 1])
fsarea = np.abs(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])) / (2. * 3456 * 5184)
fs_area.append(fsarea)

# <codecell>

# Area densities inside and outside the focus
in_area_density.append(incount * np.mean(insize) / fsarea)
out_area_density.append(outcount * np.mean(outsize) / (1. - fsarea))
                        
# Point densities inside and outside the focus
in_point_density.append(incount / fsarea)
out_point_density.append(outcount / (1. - fsarea))

# <codecell>

# Variance and CoV of outside the focus
mask = np.zeros_like(clean)
for c in cent :
    mask[int(np.round(c[0])), int(np.round(c[1]))] = 1
mask = ndimage.binary_fill_holes(ndimage.maximum_filter(mask, 35))

areadensity = filter.gaussian_filter(clean, 51)
var_area_density.append(np.var(areadensity[mask == 0]))
cov_area_density.append(np.std(areadensity[mask == 0] / np.mean(areadensity[mask == 0])))

pointdensity = np.zeros_like(clean)
for n in nuclei :
    pointdensity[int(np.round(n.centroid[0])), int(np.round(n.centroid[1]))] = 1
pointdensity = filter.gaussian_filter(pointdensity, 51)

var_point_density.append(np.var(pointdensity[mask == 0]))
cov_point_density.append(np.std(pointdensity[mask == 0] / np.mean(pointdensity[mask == 0])))

# <codecell>

# Pairwise distances between nuclei outside the focal area 
distances = distance.squareform(distance.pdist([n.centroid for n in outer_nuclei]))

for i in range(len(distances)) :
    distances[i, i] = np.inf
    
pairwise.append(np.min(distances, axis=0).mean())

# <codecell>

print in_count
print out_count
print in_size
print out_size
print in_area_density
print out_area_density
print in_point_density
print out_point_density
print var_point_density
print cov_point_density
print var_area_density
print cov_area_density
print fs_area

