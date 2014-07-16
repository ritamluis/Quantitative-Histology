# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# Files
import os

# Basic
import numpy as np

# Image Processing
import skimage as ski
from skimage import io, feature, morphology, filter, exposure, color, transform
import scipy.signal as sig

# Stats
import scipy.stats as st

# Nonlinear Fitting
import lmfit as lm

# Visualisation
import matplotlib
import matplotlib.pyplot as plt
import seaborn
#matplotlib.rcParams['savefig.dpi'] = 3. * matplotlib.rcParams['savefig.dpi']

# <codecell>


# <codecell>

def normalise(im) :
    return (im - im.min()).astype(float) / (im - im.min()).astype(float).max()

# <codecell>

# Read files
files = []

directory = "data/"
for i in os.listdir(directory) :
    if i.endswith(".jpg") :
        if not i.endswith("processed.jpg") :
            files.append(directory + i)
            
print files

# <codecell>


# <codecell>

A = io.imread(files[0])
As = transform.rescale(A, 0.25)
io.imshow(A)
plt.grid(False)

# <codecell>

#B = exposure.adjust_sigmoid(A, gain=12)
Bs = exposure.adjust_sigmoid(ski.img_as_float(As), gain=12)
#io.imshow(B - exposure.adjust_sigmoid(ski.img_as_float(A), gain=12))

# <codecell>

#C = color.rgb2xyz(B)[:, :, 1]
Cs = color.rgb2xyz(Bs)[:, :, 1]
io.imshow(Cs)
plt.grid(0)

# <codecell>

#D = filter.threshold_adaptive(C, 301)
Ds = filter.threshold_adaptive(Cs, 75)
io.imshow(Ds)
plt.grid(0)

# <codecell>

#E = morphology.remove_small_objects(~morphology.remove_small_objects(~D, 100), 100)
Es = morphology.remove_small_objects(~morphology.remove_small_objects(~Ds, 10), 10)
io.imshow(Es)
plt.grid(False)

# <codecell>

X = Es.copy()#transform.rescale(E, 0.25)

scales = 1. / np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20]) # 2 to 20 pixels
orientations = np.linspace(0, np.pi * 17./18., 18) # 0 to 180 degrees in 10 degree increments

# Results array
gabor = np.zeros((len(orientations), len(scales)))

# Perform Gabor filtering
for i, iv in enumerate(orientations) :
    for j, jv in enumerate(scales) :
        Y, Z = filter.gabor_filter(ski.img_as_float(X), jv, iv)
        gabor[i, j] = np.sqrt(np.sum(np.abs(Y)**2) + np.sum(np.abs(Z)**2)) # Return energy
        
    print i

# <codecell>

px = Y.ravel()
py = 1 - px

# <codecell>


# <codecell>



io.imshow(gabor)
plt.grid(False)

# <codecell>

plt.plot(gabor[:,5])

# <codecell>

yy = y[:, 8]
zz = z[:, 6]
params = lm.Parameters()
params.add("amp", value = yy.max() - yy.min(), min = 0, max = yy.max())
params.add("intercept", value = yy.min(), min = 0, max = yy.max())
params.add("std", value = 0.2, min = 0)
params.add("mean", value = np.where(yy == yy.max())[0][0], min = 0)

result = lm.minimize(gaussian, params, args = (orientations, yy))

# <codecell>

lm.report_fit(params)
yyy = params["amp"].value * np.exp( (-(orientations - params["mean"].value)**2) / (2*(params["std"].value**2)) ) + params["intercept"].value

# <codecell>

print (yy.max() - yy.min()) / yy.mean(), (zz.max() - zz.min()) / zz.mean()

# <codecell>

%%timeit 
xx = np.zeros((500,500))
Y, Z = filter.gabor_filter(xx, 0.05, 0)

