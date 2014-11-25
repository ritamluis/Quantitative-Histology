# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os

import numpy as np

import skimage as ski
from skimage import io, filter, color, exposure, morphology, feature, draw, measure, transform

import matplotlib
%matplotlib inline
matplotlib.rcParams["figure.figsize"] = (16, 10)

# <codecell>

l = os.listdir("../data")
print l[6]

qstain = np.array([[.26451728, .5205347, .81183386], [.9199094, .29797825, .25489032], [.28947765, .80015373, .5253158]])

# <codecell>

A = io.imread("../data/" + l[5])
io.imshow(A)

# <codecell>

#B = exposure.adjust_sigmoid(filter.gaussian_filter(A[:, :, 0], 19), cutoff=.45, gain=15)



B = exposure.adjust_sigmoid(
    filter.gaussian_filter(
        exposure.rescale_intensity(
            color.separate_stains(
                A, 
                np.linalg.inv(qstain)), 
            out_range=(0, 1))[:, :, 1], 
        29), 
    cutoff=.35, gain=20)
io.imshow(B)

# <codecell>

#b = morphology.remove_small_objects(filter.threshold_adaptive(filter.gaussian_filter(exposure.adjust_sigmoid(A[:,:,1]), 31), 501, offset=-0.05), 2000)

C = morphology.remove_small_objects(
    filter.threshold_adaptive(B, 301, offset=-0.025), 
    4000)
#io.imshow(morphology.binary_closing(np.logical_or(morphology.binary_dilation(C, morphology.disk(11)), b), morphology.disk(31)))
io.imshow(C)


#io.imshow(exposure.adjust_sigmoid(A[:, :, 1]))

# <codecell>

io.imshow(A)

# <codecell>

d = ski.img_as_float(color.separate_stains(A, np.linalg.inv(qstain)))

dd = exposure.adjust_sigmoid(filter.gaussian_filter(exposure.rescale_intensity(d[:, :, 0] + d[:, :, 2], out_range=(0, 1)), 21), gain=7, cutoff=0.6)

ddd = morphology.remove_small_objects(filter.threshold_adaptive(dd, 301, offset=-0.06), 2000)


io.imshow(morphology.binary_erosion(morphology.binary_dilation(ddd, morphology.disk(41)), morphology.disk(21)))

#io.imshow(exposure.rescale_intensity(ski.img_as_float(A[:, :, 1]) + ski.img_as_float(A[:, :, 2]), out_range=(0, 256)))

# <codecell>

e = filter.gaussian_filter(exposure.rescale_intensity(d[:,:,1], out_range=(0, 1)), 31)
#io.imshow(filter.gaussian_filter(e, 21))
#print np.mean(e)
plt.plot(exposure.histogram(e)[1], exposure.histogram(e)[0])

# <codecell>

io.imshow(filter.threshold_adaptive(e, 301, offset=0.025))

# <codecell>

io.imshow(A)

# <codecell>

from skimage import data
ihc = data.immunohistochemistry()
ihc_hdx = color.separate_stains(ihc, color.hdx_from_rgb)
ihc_rgb = color.combine_stains(ihc_hdx, color.rgb_from_hdx)

io.imshow(ihc_hdx[:,:,0])

# <codecell>

print [[0.644211, .835*0.716556, 0.266844], [0.092789, .835*0.954111, 0.283111], [0.00001, 0.00001, 0.00001]]
io.imshow(color.separate_stains(A, 
    np.linalg.inv([[0.644211, 0.716556, 0.266844], [0.092789, 0.954111, 0.283111], [0.00001, 0.00001, 0.00001]]))[:, :, 2])



# <codecell>

qstain = np.array([[.26451728, .5205347, .81183386], [.9199094, .29797825, .25489032], [.28947765, .80015373, .5253158]])
qfrompaper = np.array([[0.644211, .835*0.716556, 0.266844], [0.092789, .835*0.954111, 0.283111], [0.75919851748933231, 0.085468001712004443, 0.92121791197468583]])

aa = io.imread("/Users/qcaudron/Desktop/he copy.jpg")
ab = color.separate_stains(A,
    np.linalg.inv(qstain))

ac = color.combine_stains(ab,
    (qstain))

io.imshow(ab[1000:2000, :1000, 1])
#print np.max(ab[:,:,2])

# <codecell>

np.sum(q**2, axis=1)

# <codecell>


# <codecell>

print (color.rgb_from_hed)
print color.hed_from_rgb

