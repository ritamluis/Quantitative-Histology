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
print l
qstain = np.array([[.26451728, .5205347, .81183386], [.9199094, .29797825, .25489032], [.28947765, .80015373, .5253158]])

# <codecell>

# Original image
A = transform.rescale(io.imread("../data/" + l[0]), 0.25)
io.imshow(A)

# <codecell>

# Colour Deconvolution
deconv = ski.img_as_float(color.separate_stains(A, np.linalg.inv(qstain)))

"""
io.imshow(exposure.equalize_adapthist(
          exposure.adjust_sigmoid(
            filter.gaussian_filter(
                exposure.rescale_intensity(
                    1-A[:, :, 2], 
                    out_range=(0, 1)), 
                4), 
            gain=7, cutoff=0.7)))
"""

#io.imshow(exposure.adjust_sigmoid(filter.gaussian_filter(exposure.equalize_adapthist(exposure.rescale_intensity(1 - A[:, :, 2], out_range=(0, 1)), ntiles_y=1), 3), cutoff=0.55))



# <codecell>

# Blood

"""
blood = \
morphology.remove_small_objects(
    morphology.binary_closing(
        morphology.remove_small_objects(
            filter.threshold_adaptive(
                exposure.adjust_sigmoid(
                    filter.gaussian_filter(
                    exposure.equalize_adapthist(
                        exposure.rescale_intensity(
                            1. - A[:, :, 2], 
                            out_range=(0, 1)), 
                        ntiles_y = 1),
                        3), 
                    cutoff=0.55),
                100, offset=-0.2),
            50),
        morphology.disk(20)).astype(bool),
    100)

io.imshow(blood)
print type(blood)
"""

# <codecell>

# Veins

veins = \
morphology.remove_small_objects(
    binary_fill_holes(
        morphology.binary_closing(
            morphology.remove_small_objects(      
                filter.threshold_adaptive(
                    filter.gaussian_filter(      
                        maximum_filter(
                            deconv[:, :, 2] / deconv[:, :, 0],
                            3),
                        7),
                    170, offset=-0.13),
                60),
            morphology.disk(29))
        ),
    250)

io.imshow(veins)

# <codecell>

io.imshow( \
          filter.gaussian_filter(      
	                        maximum_filter(
                                exposure.rescale_intensity(
	                            deconv[:, :, 2] / deconv[:, :, 0], out_range=(0,1)),
	                            2),
	                        11))

# <codecell>

# Inflammation

inflammation = \
morphology.binary_closing(
    morphology.remove_small_objects(
        filter.threshold_adaptive(
            exposure.adjust_sigmoid(
                filter.gaussian_filter(
                    exposure.equalize_adapthist(
                        exposure.rescale_intensity(
                            deconv[:, :, 1],
                            out_range = (0, 1)),
                        ntiles_y = 3),
                    5),
                cutoff = 0.6),
            50, offset = -0.1),
        150),
    morphology.disk(19))



"""
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
"""
io.imshow(inflammation)
#dec = exposure.adjust_sigmoid(filter.gaussian_filter(exposure.equalize_adapthist(exposure.rescale_intensity(deconv[:, :, 1], out_range=(0, 1)), ntiles_y=3), 5), cutoff=0.6)
#io.imshow(filter.threshold_adaptive(dec, 50, offset=-0.1))
#io.imshow(filter.threshold_adaptive(exposure.adjust_sigmoid(filter.gaussian_filter(exposure.rescale_intensity(deconv[:, :, 1], out_range=(0, 1)), 19), gain=15, cutoff=0.4), 401))

# <codecell>

# Labelled
total = np.zeros_like(A)
#total[:, :, 0] = blood
total[:, :, 1] = veins
total[:, :, 2] = inflammation

io.imshow(total)

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

#full = (np.logical_or(np.logical_or(blood, veins), inflammation))

#io.imshow(binary_fill_holes(minimum_filter(maximum_filter(full, size=51, footprint=morphology.diamond(31)), size=51)))

# <codecell>

#io.imshow(ski.img_as_float(A)*0.5 + full[:, :, newaxis]*0.5)
#io.imshow(morphology.binary_closing(full, morphology.disk(15)))

# <codecell>

io.imsave("__1.gif", A)
io.imsave("__2.gif", total)
#io.imsave("__2.gif", (morphology.binary_closing(full, morphology.disk(15))).astype(float))

# <codecell>

os.system("gifsicle --delay=80 --loop __*.gif > test.gif")

# <codecell>

print B.shape
print transform.rescale(B, 0.5).shape

# <codecell>


