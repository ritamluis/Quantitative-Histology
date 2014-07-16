# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os

import numpy as np

from scipy.ndimage import maximum_filter, minimum_filter, binary_fill_holes

import skimage as ski
from skimage import io, filter, color, exposure, morphology, feature, draw, measure, transform

figsize(16, 10)

# <codecell>

l = os.listdir("../data")
print l[6]

qstain = np.array([[.26451728, .5205347, .81183386], [.9199094, .29797825, .25489032], [.28947765, .80015373, .5253158]])

# <codecell>

# Original image
A = io.imread("../data/" + l[2])
io.imshow(A)

# <codecell>

# Colour Deconvolution
deconv = ski.img_as_float(color.separate_stains(A, np.linalg.inv(qstain)))

io.imshow(exposure.adjust_sigmoid(
            filter.gaussian_filter(
                exposure.rescale_intensity(
                    deconv[:, :, 0] + deconv[:, :, 2], 
                    out_range=(0, 1)), 
                17), 
            gain=7, cutoff=0.6))

# <codecell>

# Blood

blood = \
morphology.remove_small_objects(
    filter.threshold_adaptive(
        exposure.adjust_sigmoid(
            filter.gaussian_filter(
                exposure.rescale_intensity(
                    deconv[:, :, 0] + deconv[:, :, 2], 
                    out_range=(0, 1)), 
                15), 
            gain=10, cutoff=0.7),
        301, offset=-0.15),
    500)

io.imshow(blood)

# <codecell>

# Veins

veins = \
binary_fill_holes(
    morphology.remove_small_objects(
        filter.threshold_adaptive(
            filter.gaussian_filter(
                exposure.adjust_sigmoid(A[:,:,1]), 
                31), 
            501, offset=-0.07), 
        3000) #### SRSLY ?
    )

io.imshow(veins)

# <codecell>

# Inflammation

inflammation = \
morphology.remove_small_objects(
    filter.threshold_adaptive(
        exposure.adjust_sigmoid(
            filter.gaussian_filter(
                exposure.rescale_intensity(
                    deconv[:, :, 1],
                    out_range = (0, 1)),
                25),
            gain = 12, cutoff = 0.32),
        501, offset = -0.1),
    4000)


io.imshow(inflammation)
#io.imshow(filter.threshold_adaptive(exposure.adjust_sigmoid(filter.gaussian_filter(exposure.rescale_intensity(deconv[:, :, 1], out_range=(0, 1)), 19), gain=15, cutoff=0.4), 401))

# <codecell>

# Labelled
io.imshow(veins + 2 * blood + 3 * inflammation)

# <codecell>

# All regions of interest

#%timeit

"""
all = \
minimum_filter(
    maximum_filter(
        np.logical_or(np.logical_or(blood, veins), inflammation).astype(float),
        size = 101),
    size = 51)
"""



"""
morphology.binary_erosion(
    morphology.binary_dilation(
        np.logical_or(np.logical_or(blood, veins), inflammation),
        morphology.disk(31)),
    morphology.disk(21))
"""

full = (np.logical_or(np.logical_or(blood, veins), inflammation))

#io.imshow(binary_fill_holes(minimum_filter(maximum_filter(full, size=51, footprint=morphology.diamond(31)), size=51)))

# <codecell>

io.imshow(ski.img_as_float(A)*0.5 + full[:, :, newaxis]*0.5)

# <codecell>

io.imsave("__1.gif", A)
io.imsave("__2.gif", ski.img_as_float(full))

# <codecell>

os.system("gifsicle --delay=80 --loop __*.gif > test.gif")

# <codecell>


