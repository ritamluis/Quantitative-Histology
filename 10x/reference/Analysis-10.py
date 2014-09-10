# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os

import numpy as np

import skimage as ski
from skimage import io, filter, color, exposure, morphology, feature, draw, measure, transform

from sklearn.linear_model import BayesianRidge
from sklearn.mixture import GMM

from scipy import spatial, ndimage, signal, stats

import matplotlib.pyplot as plt
import seaborn
figsize(16, 10)

#%matplotlib inline


# WANT
# - Density outside inflammatory zone
# - Area, width, length of focus
# - Area of inflammation / area of vein

# <codecell>

# Load filenames
dirq = "/Users/qcaudron/repositories/Quantitative-Histology/10x/data/"
files = []
for i in os.listdir(dirq) :
    if i.endswith(".jpg") :
        files.append(dirq + i)
        

offset = 1000
lengthscale = 301



image = exposure.equalize_adapthist(io.imread("/Users/qcaudron/repositories/Quantitative-Histology/10x/data/Sheep24-10x-19.jpg"))#files[4]))
#image = io.imread(files[4])
io.imshow(image)
plt.grid(False)

# <codecell>

#binary = np.logical_or(filter.threshold_adaptive(exposure.adjust_sigmoid(image[:, :, 0], cutoff=0.4, gain=20), lengthscale), \
#                       filter.threshold_adaptive(exposure.adjust_sigmoid(image[:, :, 2], cutoff=0.5, gain=20), lengthscale))

binary = filter.threshold_adaptive(exposure.adjust_sigmoid(image[:, :, 0], cutoff=0.4, gain = 30), 301).astype(bool)
clean = morphology.binary_closing(binary, morphology.disk(3)).astype(bool)
clean = morphology.remove_small_objects(clean, 200)
clean = morphology.remove_small_objects( (1-clean).astype(bool), 200)



io.imshow(clean)
plt.grid(False)


# <codecell>

(xdim, ydim, _) = image.shape
xdim /= 10
ydim /= 10

local_density = filter.gaussian_filter(clean, 61)

fig, ax = plt.subplots()

io.imshow(local_density)

#for n, c in enumerate(measure.find_contours(local_density, local_density.max()/2)) :
#    ax.plot(c[:, 1], c[:, 0], linewidth=2)

plt.grid(False)

# <codecell>


ent = filter.gaussian_filter(filter.rank.entropy(local_density, morphology.disk(3)), 75)

ent -= ent.min()
ent /= ent.max()

local_density -= local_density.min()
local_density /= local_density.max()

io.imshow(2 * ent + local_density)
plt.grid(False)

# <codecell>

info = ent * (1 + local_density)

io.imshow(info)
plt.grid(False)

# <codecell>

bw = (info) > filter.threshold_otsu(info)
#bw = morphology.dilation(bw, morphology.disk(21))


fig, ax = plt.subplots()

io.imshow(image)
C = measure.find_contours(bw, 0.5)
centroid = []
vals = []
for c in C :
#    ax.plot(c[:, 1], c[:, 0], lw=5)
    centroid.append(np.linalg.norm([c[:, 1].mean() - bw.shape[1] / 2, c[:, 0].mean() - bw.shape[0] / 2]))
    vals.append(local_density.T[c.astype(int)].sum())
cent = C[np.argmin(centroid / np.array(vals))]
    
ax.plot(cent[:, 1], cent[:, 0], lw=5, c="k", alpha = 0.7)

plt.grid(False)

# <codecell>

plt.scatter(centroid, centroid / np.array(vals))
#ax = plt.gca()
#ax.set_yscale("log")
#ax.set_ylim([1e-9, 1e-5])

print 1. / np.array(vals)

# <codecell>

# DO NOT RUN ME

for f in files :
    image = exposure.equalize_adapthist(io.imread(f))
    binary = filter.threshold_adaptive(exposure.adjust_sigmoid(image[:, :, 0], cutoff=0.4, gain = 30), 301).astype(bool)
    clean = morphology.binary_closing(binary, morphology.disk(3)).astype(bool)
    clean = morphology.remove_small_objects(clean, 200)
    clean = morphology.remove_small_objects( (1-clean).astype(bool), 200)
    local_density = filter.gaussian_filter(clean, 61)
    ent = filter.gaussian_filter(filter.rank.entropy(local_density, morphology.disk(3)), 101)
    ent -= ent.min()
    ent /= ent.max()
    local_density -= local_density.min()
    local_density /= local_density.max()
    info = ent * (1 + local_density)
    bw = (info) > filter.threshold_otsu(info)


    fig, ax = plt.subplots()
    plt.imshow(image)
    C = measure.find_contours(bw, 0.5)
    centroid = []
    for c in C :
        centroid.append(np.linalg.norm([c[:, 1].mean() - bw.shape[1] / 2, c[:, 0].mean() - bw.shape[0] / 2]))
    cent = C[np.argmin(centroid)]

    ax.scatter(cent[:, 1], cent[:, 0], lw=3)

    plt.grid(False)
    plt.savefig(f.split("Sheep")[1])
    print ("Done " + f.split("Sheep")[1])

# <codecell>

f.split("Sheep")

# <codecell>


# <codecell>

# Segmentation

distance = ndimage.distance_transform_edt(clean)
peaks = feature.peak_local_max(distance, indices=False, labels = clean)
markers = ndimage.label(peaks)[0]
labels = morphology.watershed(-distance, markers, mask = clean)

io.imshow(labels, interpolation = "nearest", cmap = plt.cm.spectral)

# <codecell>

io.imshow(filter.gaussian_filter(image, 31))

# <codecell>

# Local density

B = image[:, :, 2]
#B = filter.gaussian_filter(np.abs(B - B.mean()), 31)
#io.imshow(B > 1.1*filter.threshold_otsu(B))
io.imshow(color.rgb2xyz(image)[:, :, 2])


#local_density = filter.gaussian_filter(image[:,:,0], 41)# - filter.gaussian_filter(clean, 201)
#io.imshow(local_density > 1.2 * filter.threshold_otsu(local_density, nbins = 1000))
#io.imshow(local_density - local_density.mean() > 0)
"""
#io.imshow(filter.gaussian_filter(clean, 51))
Q = (signal.convolve2d(ski.img_as_float(clean / clean.max()), ski.img_as_float(morphology.disk(201)), "valid"))# - filter.gaussian_filter(clean, 251))
Q -= Q.min()
Q /= Q.max()
io.imshow(Q)
"""

# <codecell>

B = signal.fftconvolve(ski.img_as_float(A), morphology.disk(31))
B -= B.min()
B /= B.max()
isignal.fftconvolve

# <codecell>

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

B = np.zeros((100,100))
B[30:70, 30:70] = 1
using_c2d = signal.convolve2d(B, B, mode = "same")
using_fft = signal.fftconvolve(B, B, mode = "same")

plt.subplot(131)
plt.imshow(B)
plt.grid(False)

plt.subplot(132)
plt.imshow(using_c2d)
plt.grid(False)

plt.subplot(133)
plt.imshow(using_fft)
plt.grid(False)

print ((using_fft - using_c2d) < 1e-10).all()

# <codecell>

plt.imshow(using_fft - using_c2d)

# <codecell>

def fftconvolve2d(x, y, mode="full"):
    """
    x and y must be real 2-d numpy arrays.

    mode must be "full" or "valid".
    """
    x_shape = np.array(x.shape)
    y_shape = np.array(y.shape)
    z_shape = x_shape + y_shape - 1
    z = ifft2(fft2(x, z_shape) * fft2(y, z_shape)).real

    if mode != "full":
        # To compute a valid shape, either np.all(x_shape >= y_shape) or
        # np.all(y_shape >= x_shape).
        valid_shape = x_shape - y_shape + 1
        if np.any(valid_shape < 1):
            valid_shape = y_shape - x_shape + 1
            if np.any(valid_shape < 1):
                raise ValueError("empty result for valid shape")

        if mode == "valid" :                    
            start = (z_shape - valid_shape) // 2
            end = start + valid_shape
            z = z[start[0]:end[0], start[1]:end[1]]
            
        elif mode == "same" :
            start = (z_shape - x_shape) // 2
            z = z[start[0]:-start[0], start[1]:-start[1]]
        
        

    return z

def makeGaussian(size, fwhm = 3, center=None):
    """ Make a square gaussian kernel.
 
    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """
 
    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]
    
    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]
    
    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)

A = clean[1000:2000, 2500:]
io.imshow(fftconvolve2d(A, makeGaussian(61, fwhm=31), mode="same"))

# <codecell>


# <codecell>

for im in files :
    
    # Read image
    image = io.imread(im)
    
    # Threshold for nuclei
    binary = np.logical_or(filter.threshold_adaptive(exposure.adjust_sigmoid(image[:, :, 0], cutoff=0.4, gain=30), lengthscale), \
                      filter.threshold_adaptive(exposure.adjust_sigmoid(image[:, :, 2], gain=30), lengthscale) )
    
    # Clean up binary image
    clean = (1 - morphology.binary_closing(morphology.remove_small_objects(binary, lengthscale), morphology.disk(3))) #* \
#            color.rgb2grey(image)


    io.imshow(image)
    
    # Create smooth density surface
#    surf = filter.gaussian_filter(clean, lengthscale, mode = "constant")
    
    # Calculated baseline and remove
    """
    r = [BayesianRidge().fit(np.expand_dims(np.arange(surf.shape[i]), axis=1), np.mean(surf, axis = 1-i)) for i in [0, 1]]
    X, Y = np.meshgrid(r[1].predict(np.expand_dims(np.arange(surf.shape[1]), axis=1)), \
                       r[0].predict(np.expand_dims(np.arange(surf.shape[0]), axis=1)))
    """


"""
    r = [BayesianRidge().fit(np.vander(np.arange(surf.shape[i]), 2), np.mean(surf, axis = 1-i)) for i in [0, 1]]
    X, Y = np.meshgrid(r[1].predict(np.vander(np.arange(surf.shape[1]), 2)), \
                       r[0].predict(np.vander(np.arange(surf.shape[0]), 2)))
    
    baseline = X * Y
#    hist, bins = np.histogram(surf.ravel(), 100)
#    baseline = bins[np.argmax(hist)+1]# + np.std(surf)
    zeroed = surf - baseline
    zeroed[zeroed < 0] = 0
    
    # Find peaks in image centre
    peaks = feature.peak_local_max(zeroed[offset:-offset, offset:-offset])
    dist = np.min(spatial.distance.pdist(peaks))/2. if len(peaks) > 1 else min(zeroed.shape)/2.
    
    # Compute volume integrals and maxima per peak
    integrals = []
    maxima = []
    for i in peaks :
        mask = np.zeros_like(zeroed)
        mask[draw.circle(i[0] + offset, i[1] + offset, dist, shape = zeroed.shape)] = 1
        integrals.append(np.sum(zeroed * mask) / np.sum(mask))
        maxima.append((zeroed[i[0] + offset, i[1] + offset] + baseline) / baseline)
"""    

#    io.imshow(surf)
#    plt.grid(0)

# <codecell>

"""
#image_max = ndimage.maximum_filter(dist, size=11, mode='constant')


#dist = ndimage.distance_transform_edt(clean)
#crit = np.zeros_like(clean)
#coordinates = feature.peak_local_max(dist, min_distance=5)

s2 = filter.gaussian_filter(clean, 21)
dist2 = ndimage.distance_transform_edt(s2 > filter.threshold_otsu(s2))

coordinates2 = feature.peak_local_max(dist2, min_distance=5)

#print len(coordinates), len(coordinates2)
#scales = [11, 21, 31, 41, 51]

spatial.voronoi_plot_2d(spatial.Voronoi(coordinates2))
#io.imshow(filter.gaussian_filter(clean, 11))
"""


nuclei = feature.peak_local_max(ndimage.distance_transform_edt(clean), min_distance=11)
tessellation 

scales = [11, 21, 31, 41, 51]
for i in scales :

# <codecell>

x = []
for i in range(3, 51, 2) :
    s = filter.gaussian_filter(clean, i)
    x.append(len(feature.peak_local_max(ndimage.distance_transform_edt(s > filter.threshold_otsu(s)), min_distance=5)))
    
plt.plot(range(3, 51, 2), x)

# <codecell>

Q = [clean]
kernelsize = 51
downsamplesize = 3
for i in range(1, 6) :
    im = filter.gaussian_filter(morphology.binary_opening(Q[i-1], morphology.disk(5)), kernelsize)
    im2 = zeros((im.shape[0]/downsamplesize, im.shape[1]/downsamplesize))
    for x in range(im2.shape[0]) :
        for y in range(im2.shape[1]) :
            im2[x, y] = np.mean(im[downsamplesize*x:downsamplesize*(x+1), downsamplesize*y:downsamplesize*(y+1)])
    im2 -= im2.min()
    im2 /= im2.max()
    Q.append(np.round(im2*1.))

# <codecell>

q = []
for i in Q :
    q.append(np.sum(i) / (i.shape[0] * i.shape[1]))
plt.plot(q)

# <codecell>

io.imshow(Q[2])

# <codecell>

#io.imshow(morphology.binary_erosion(clean, np.ones((21,21))))
#io.imshow(morphology.binary_dilation(clean, morphology.disk(25)))
io.imshow(morphology.binary_dilation(filter.rank.mean(clean, morphology.disk(51)), morphology.disk(51)))

# <codecell>

blah = filter.gaussian_filter(clean, 401)
plt.hist(blah.ravel(), 256);

# <codecell>

io.imshow(blah)# > np.median(blah.ravel()))

# <codecell>

#r = [BayesianRidge().fit(np.expand_dims(np.arange(surf.shape[i]), axis=1), np.mean(surf, axis = 1-i)) for i in [0, 1]]
X, Y = np.meshgrid(r[0].predict(np.expand_dims(np.arange(surf.shape[0]), axis=1)), \
                   r[1].predict(np.expand_dims(np.arange(surf.shape[1]), axis=1)))

zeroed.shape

# <codecell>

c = clean - np.min(clean)
c /= c.max()
c = c.astype(bool)
io.imsave("/Users/qcaudron/Desktop/charo/2_smoothed.jpg", ski.img_as_uint(surf))

# <codecell>

z1 = np.mean(surf, axis=0)
z2 = np.mean(surf, axis=1)

#for i in range(surf.shape[1]) : 
#    plt.plot(surf[:, i], "k")
#plt.plot(z2)
r = [BayesianRidge().fit(np.vander(np.arange(surf.shape[i]), 2), np.mean(surf, axis = 1-i)) for i in [0, 1]]
r1 = BayesianRidge().fit(np.arange(len(z1)).reshape(len(z1),1), z1)
r2 = BayesianRidge().fit(np.arange(len(z2[500:-500])).reshape(len(z2[500:-500]),1), z2[500:-500])

#plt.plot(r1.predict(np.arange(len(z1)).reshape(len(z1),1)), linewidth=5)
plt.plot(r2.predict(np.arange(len(z2)).reshape(len(z2),1)), linewidth=5)
plt.plot(z2, linewidth=5)
#plt.axhline(b[np.argmax(h)], c="r", linewidth=3)
#plt.plot(r[0].predict(np.vander(np.arange(surf.shape[0]), 2)), linewidth=3)

#plt.plot(r[0].predict(np.arange(len(z1)).reshape(len(z1),1)), linewidth=3)
#plt.plot(r[0].predict(np.expand_dims(np.arange(surf.shape[0]), axis=1)), linewidth=5)
#plt.axhline(np.mean(z1 / r1.predict(np.arange(len(z1)).reshape(len(z1),1))))

# <codecell>

lz = np.log(z2)
r3 = BayesianRidge().fit(np.arange(len(lz[500:-500])).reshape(len(lz[500:-500]),1), lz[500:-500])

plt.plot(np.exp(lz))
plt.plot(np.exp(r3.predict(np.arange(len(lz)).reshape(len(lz),1))))

# <codecell>

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x, y = np.meshgrid(range(surf.shape[0]), range(surf.shape[1]))

x.shape
ax.plot_surface(y, x, zeroed.T, rstride = 50, cstride = 50, antialiased = False, cmap = cm.coolwarm, linewidth=0)

# <codecell>

np.min(baseline)

# <codecell>

x2, y2 = draw.circle(peaks[1][0]+1000, peaks[1][1]+1000, dist*1., shape = (zeroed.shape))
mask2 = zeros_like(zeroed)
mask2[x, y] = 1
np.sum(zeroed * mask - zeroed * mask2)

np.sum(zeroed * mask2)

