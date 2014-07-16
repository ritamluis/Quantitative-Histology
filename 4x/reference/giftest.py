import os

import numpy as np

from scipy.ndimage import maximum_filter, minimum_filter, binary_fill_holes

import skimage as ski
from skimage import io, filter, color, exposure, morphology, feature, draw, measure, transform

#figsize(16, 10)

l = os.listdir("../data")
for f in l :
	if not f.endswith(".jpg") :
		l.remove(f)


qstain = np.array([[.26451728, .5205347, .81183386], [.9199094, .29797825, .25489032], [.28947765, .80015373, .5253158]])




for im in l :

	print im

	A = transform.rescale(io.imread("../data/" + im), 0.25)

	deconv = ski.img_as_float(color.separate_stains(A, np.linalg.inv(qstain)))






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





	
	inflammation = \
	maximum_filter(
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
		    250),
		29)







	# Labelled
	total = np.zeros_like(A)
	#total[:, :, 0] = blood
	total[:, :, 1] = veins
	total[:, :, 2] = inflammation







	io.imsave("__1.gif", A)
	io.imsave("__2.gif", total)

	os.system("gifsicle --delay=80 --loop __*.gif > %s.gif" % im)

	