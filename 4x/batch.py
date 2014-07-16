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
from scipy.ndimage import maximum_filter, minimum_filter, binary_fill_holes

# Stats
import scipy.stats as st

# Visualisation
#import matplotlib
#import matplotlib.pyplot as plt
#import seaborn




# Results : sheep ID, entropy, entropic variance, lacunarity, directionality
# IDs
name = []
imageid = []

# Entropy
ent = []
entstd = []

# Lacunarity
lac = []
nonlac = []

# Gabor filtering
direc = []
scale = []

# Inflammatory foci count and size
inflcount = []
inflsize = []
inflstd = []

# Tissue to sinusoid ratio
ratio = []








# Image files to process

# If we're running on all files
if len(sys.argv) == 1 :

    # Read files
    files = []
    directory = "data/"
    for i in os.listdir(directory) :
        if i.endswith(".jpg") : # if it's a jpg
            if not i.endswith("_processed.jpg") : # and isn't a processed image
                files.append(directory + i) # then add it to the list to be processed


# Otherwise, we ran launcher.py and need to do select images only
else :

    F = open("filelist%i.p" % int(sys.argv[1]), "r")
    files = pickle.load(F)
    F.close()











# Iterate over files
for f in files :

    print "Thread %s. Opening file : %s" % (sys.argv[1], f)


    # If it's not been processed before
    if not os.path.exists(f + "_processed.jpg") :

        ### PROCESSING
        
        # Read the image
        A = io.imread(f)
        
        # Constrast enhancement
        print "Thread %s. Sigmoid transform for contrast." % sys.argv[1]
        B = exposure.adjust_sigmoid(A, gain=12)
        
        # Extract luminosity
        print "Thread %s. Generating luminosity." % sys.argv[1]
        C = color.rgb2xyz(B)[:, :, 1]
        
        # Apply adaptive thresholding
        print "Thread %s. Performing adaptive thresholding." % sys.argv[1]
        D = filter.threshold_adaptive(C, 301)
        
        # Clean
        print "Thread %s. Cleaning image." % sys.argv[1]
        E = morphology.remove_small_objects(~morphology.remove_small_objects(~D, 100), 100)

        # Save to disk
        io.imsave(f + "_processed.jpg", ski.img_as_float(E))
        
        # Downsample for Gabor filtering
        print "Thread %s. Downsampling image." % sys.argv[1]
        Es = ski.img_as_float(transform.rescale(E, 0.25))


    else :

        # Otherwise, we've processed it before, so read it in for speed
        A = io.imread(f)
        E = ski.img_as_float(io.imread(f + "_processed.jpg"))
        Es = ski.img_as_float(transform.rescale(E, 0.25))
    print "Thread %s. Image already processed, downsampled." % sys.argv[1]
        

    ## ID
    name.append(int(f.split("data/Sheep")[1].split("-")[0].split(".jpg")[0]))
    imageid.append(int(f.split("data/Sheep")[1].split("-")[2].split(".jpg")[0]))














    

    

    ### GABOR FILTERING

    # Define scales and orientations to compute over
    pixelscales = np.arange(15, 55, 4)
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

    # Results
    scale.append(pixelscales[optimalscale])
    direc.append(directionality)



    












    ### LACUNARITY AND ENTROPY

    # Define scales over which to compute lacunarity and entropy
    
    scaleent = []
    scaleentstd = []
    #scalelac = []
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
    #px = np.reshape(Y[int(i) : int(E.shape[0]-i), int(i) : int(E.shape[1]-i)] / circlesize, (E.shape[0]-2*i)*(E.shape[1]-2*i),)
    py = 1. - px
    idx = np.logical_and(px > 1. / circlesize, px < 1.)
    entropy = - ( np.mean(px[idx] * np.log(px[idx])) + np.mean(py[idx] * np.log(py[idx])) )
    entropystd = np.std(px[idx] * np.log(px[idx]) + py[idx] * np.log(py[idx]))

    # Compute normalised lacunarity
    lx = np.var(Y) / (np.mean(Y) ** 2) + 1
    ly = np.var(1. - Y) / (np.mean(1. - Y) ** 2) + 1
    # lacunarity = 2. - (1. / lx + 1. / ly)

    # Roy et al, J. Struct. Geol. 2010

    # Results
    roylac.append((lx - 1.) / (1./np.mean(E) - 1.))
    scaleent.append(entropy)
    scaleentstd.append(entropystd)
    # scalelac.append(lacunarity)

    print "Thread %s. Calculated entropy and lacunarity." % sys.argv[1]




    # Old code evaluates lacunarity and entropy at all scales
    # We only need it at the image's characteristic scale
    """
    # Iterate over scales
    for i, s in enumerate(lacscales) :

        # Generate a disk at this scale
        circle = morphology.disk(s)
        circlesize = circle.sum()

        # Convolve with image
        Y = sig.fftconvolve(E, circle, "valid")

        # Compute information entropy
        px = Y.ravel() / circlesize
        #px = np.reshape(Y[int(i) : int(E.shape[0]-i), int(i) : int(E.shape[1]-i)] / circlesize, (E.shape[0]-2*i)*(E.shape[1]-2*i),)
        py = 1. - px
        idx = np.logical_and(px > 1. / circlesize, px < 1.)
        entropy = - ( np.mean(px[idx] * np.log(px[idx])) + np.mean(py[idx] * np.log(py[idx])) )
        entropystd = np.std(px[idx] * np.log(px[idx]) + py[idx] * np.log(py[idx]))

        # Compute normalised lacunarity
        lx = np.var(Y) / (np.mean(Y) ** 2) + 1
        ly = np.var(1. - Y) / (np.mean(1. - Y) ** 2) + 1
        # lacunarity = 2. - (1. / lx + 1. / ly)

        # Roy et al, J. Struct. Geol. 2010

        # Results
        roylac.append((lx - 1.) / (1./np.mean(E) - 1.))
        scaleent.append(entropy)
        scaleentstd.append(entropystd)
        # scalelac.append(lacunarity)

        print "Thread %s. Calculating entropy and lacunarity. Completion : %f" % (sys.argv[1], (float(i) / len(lacscales)))
    """

    nonlac.append(roylac)
    ent.append(scaleent)
    entstd.append(scaleentstd)
    # lac.append(scalelac)




















    # Inflammatory focus count

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
        27)






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
        29)

   

    print "Thread %s. Image segmentation complete." % sys.argv[1]




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



    inflcount.append(inflammationcount)
    inflsize.append(np.mean(inflammationsize))
    inflstd.append(np.std(inflammationsize))










    # Tissue to sinusoid ratio

    tissue = np.sum(ski.img_as_bool(Es) * (1 - ski.img_as_bool(total)))
    sinusoids = np.sum(Es.shape[0] * Es.shape[1] - np.sum(ski.img_as_bool(total)))
    ratio.append(tissue.astype(float) / sinusoids)








# Pickle results

results = {}
results["ratio"] = ratio
results["inflcount"] = inflcount
results["inflsize"] = inflsize
results["inflstd"] = inflstd
results["entropy"] = ent
results["entstd"] = entstd
results["directionality"] = direc
results["scale"] = scale
results["name"] = name
results["ID"] = imageid
results["roylac"] = nonlac

F = open("results/thread_%s" % sys.argv[1], "w")
pickle.dump(results, F)
F.close()















