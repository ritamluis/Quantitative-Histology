# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# Files
import os
import pickle
import sys

# Basic
import numpy as np

# Image Processing
import skimage as ski
from skimage import io, feature, morphology, filter, exposure, color, transform, measure
import scipy.signal as sig
from scipy.spatial import distance
from scipy.ndimage import maximum_filter, minimum_filter, binary_fill_holes

# Stats
import scipy.stats as st

import matplotlib
%matplotlib inline
import matplotlib.pyplot as plt
matplotlib.rcParams["figure.figsize"] = (16, 10)

# <codecell>

# Read files
files = []
directory = "../data/"
for i in os.listdir(directory) :
    if i.endswith(".jpg") : # if it's a jpg
        if not i.endswith("_processed.jpg") : # and isn't a processed image
            files.append(directory + i) # then add it to the list to be processed
            
files = [files[50]]
io.imshow(io.imread(files[0]))

# <codecell>


# Iterate over files
for f in files :

    # If it's not been processed before
    if not os.path.exists(f + "_processed.jpg") :

        ### PROCESSING
        
        # Read the image
        A = io.imread(f)
        
        # Constrast enhancement
        B = exposure.adjust_sigmoid(A, gain=12)
        
        # Extract luminosity
        C = color.rgb2xyz(B)[:, :, 1]
        
        # Apply adaptive thresholding
        D = filter.threshold_adaptive(C, 301)
        
        # Clean
        E = morphology.remove_small_objects(~morphology.remove_small_objects(~D, 100), 100)

        # Save to disk
        io.imsave(f + "_processed.jpg", ski.img_as_float(E))
        
        # Downsample for Gabor filtering
        Es = ski.img_as_float(transform.rescale(E, 0.25))


    else :

        # Otherwise, we've processed it before, so read it in for speed
        A = io.imread(f)
        E = ski.img_as_float(io.imread(f + "_processed.jpg"))
        Es = ski.img_as_float(transform.rescale(E, 0.25))
        

# <codecell>

A = io.imread(files[0])[:2000, 1000:-1000, :]
A2 = filter.gaussian_filter(A, 5)
B1 = exposure.adjust_sigmoid(A, gain=12)
B2 = exposure.adjust_sigmoid(A2, gain=12)
C1 = color.rgb2xyz(B1)[:, :, 1]
C2 = color.rgb2xyz(B2)[:, :, 1]
D1 = filter.threshold_adaptive(C1, 301)
D2 = filter.threshold_adaptive(C2, 301)
E1 = morphology.remove_small_objects(~morphology.remove_small_objects(~D1, 100), 100)
E2 = morphology.remove_small_objects(~morphology.remove_small_objects(~D2, 100), 100)

d1 = filter.threshold_adaptive(C1, 301, offset=-0.01)
d2 = filter.threshold_adaptive(C2, 301, offset=-0.01)
e1 = morphology.remove_small_objects(~morphology.remove_small_objects(~d1, 100), 100)
e2 = morphology.remove_small_objects(~morphology.remove_small_objects(~d2, 100), 100)

print np.abs(E1 - e1).sum()
print np.abs(E2 - e2).sum()

# <codecell>

(A2.shape[0] * A2.shape[1])

# <codecell>

    pixelscales = np.arange(15, 55, 2)
    gaborscales = 4. / pixelscales # 2 to 20 pixels

    orientations = np.linspace(0, np.pi * 11./12., 12) # 0 to 180 degrees in 15 degree increments

    # Results array
    gabor = np.zeros((len(orientations), len(gaborscales)))

    # Perform Gabor filtering
    for i, iv in enumerate(orientations) :
        for j, jv in enumerate(gaborscales) :
            gaborReal, gaborImag = filter.gabor_filter(Es, jv, iv)
            gabor[i, j] = np.sqrt(np.sum(np.abs(gaborReal) ** 2) + np.sum(np.abs(gaborImag) ** 2)) # Return energy
        print "Thread %s. Gabor filtering. Completion : %f" % (sys.argv[1], (i / float(len(orientations))))

    # Determine orientation-independent scale which fits best
    optimalscale = np.argmax(np.sum(gabor, axis = 0))

    # At this scale, calculate directionality coefficient
    g = gabor[:, optimalscale]
    directionality = (g.max() - g.min()) / g.max()

# <codecell>

    scaleent = []
    scaleentstd = []
    roylac = []


    # Characteristic scale
    s = pixelscales[optimalscale]

    # Generate a disk at this scale
    circle = morphology.disk(s)
    circlesize = circle.sum()

    # Convolve with image
    Y = sig.fftconvolve(E, circle, "valid")

    # Compute information entropy
    px = Y.ravel() / circlesize
    py = 1. - px
    idx = np.logical_and(px > 1. / circlesize, px < 1.)
    entropy = - ( np.mean(px[idx] * np.log(px[idx])) + np.mean(py[idx] * np.log(py[idx])) )
    entropystd = np.std(px[idx] * np.log(px[idx]) + py[idx] * np.log(py[idx]))

    # Compute normalised lacunarity
    lx = np.var(Y) / (np.mean(Y) ** 2) + 1
    ly = np.var(1. - Y) / (np.mean(1. - Y) ** 2) + 1

    # Roy et al, J. Struct. Geol. 2010

    # Results
    roylac.append((lx - 1.) / (1./np.mean(E) - 1.))
    scaleent.append(entropy)
    scaleentstd.append(entropystd)

# <codecell>

    qstain = np.array([[.26451728, .5205347, .81183386], [.9199094, .29797825, .25489032], [.28947765, .80015373, .5253158]])

    deconv = ski.img_as_float(color.separate_stains(transform.rescale(A, 0.25), np.linalg.inv(qstain)))



    subveins1 = \
    morphology.remove_small_objects(
        filter.threshold_adaptive(
            filter.gaussian_filter(
                deconv[:, :, 2] / deconv[:, :, 0],
                11),
            250, offset = -0.13),
        60)

    subveins2 = \
    morphology.remove_small_objects(
        filter.threshold_adaptive(
            filter.gaussian_filter(
                    maximum_filter(
                    deconv[:, :, 2] / deconv[:, :, 0],
                    5),
                11),
            250, offset = -0.13),
        60)

    veins = \
    maximum_filter(
        morphology.remove_small_objects(
            binary_fill_holes(
                morphology.binary_closing(
                    np.logical_or(subveins1, subveins2),
                    morphology.disk(25)),
                ),
            250),
        55)






    rawinflammation = \
    morphology.remove_small_objects(
    filter.threshold_adaptive(
        exposure.adjust_sigmoid(
            filter.gaussian_filter(
                exposure.equalize_adapthist(
                    exposure.rescale_intensity(
                        deconv[:, :, 1],
                        out_range = (0, 1)),
                    ntiles_y = 1),
                    5),
            cutoff = 0.6),
            75, offset = -0.12),
        250)

    
    inflammation = \
    maximum_filter(
        rawinflammation,
        55)

# <codecell>

    total = veins + inflammation
    coloured = np.zeros_like(deconv)
    coloured[:, :, 1] = veins
    coloured[:, :, 2] = inflammation


    labelled, regions = measure.label(total, return_num = True)
    

    inflammationcount = 0
    inflammationsize = []



    for b in range(1, regions) :

        if (inflammation * (labelled == b)).sum() / ((veins * (labelled == b)).sum() + 1) > 0.5 :
            inflammationcount += 1
            inflammationsize.append((rawinflammation * labelled == b).sum())

# <codecell>

regions

# <codecell>

io.imshow(A)

# <codecell>

io.imshow(inflammation)
plt.scatter([qq[i].centroid[1] for i in range(regions-1)], [qq[i].centroid[0] for i in range(regions-1)])

# <codecell>

pairwise.min(axis=0).mean()

