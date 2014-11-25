# Collate data from numerous batch.py runs 
# into one environment and into .csv files

import copy
import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn
mpl.rcParams['text.usetex'] = True  
from scipy.stats import mode, ks_2samp, mannwhitneyu
from sklearn.neighbors import KernelDensity



scalefactor = 20.
xdim = 1.1
ydim = 1.

colours = seaborn.color_palette("deep", 8)


# Read files in results directory
files = []
newfiles = []
directory = "results/"
for f in os.listdir(directory) :
	if f.startswith("thread_") : # if it's a pickled result structure
		files.append(directory + f)
		newfiles.append(directory + "threads/" + f)

print "Found %d results files :" % len(files)
for f in files :
	print f






# Import results
results = []
for f in files :
	F = open(f)
	results.append(pickle.load(F))
	F.close()




# Process results :

nonnormlac = []
entropy = []
entstd = []
directionality = []
imageid = []
id = []
scale = []
ratio = []
inflsize = []
inflcount = []
inflstd = []
mindist = []
ifdist = []
blur = []



# Generate list of sheep IDs
uniquenames = []
for i in results :
	uniquenames.append(i["name"])

# Sort and filter
uniquenames = np.unique(np.hstack(uniquenames))






# Iterate over all unique names in order
for sheep in uniquenames :

	# Go over the results by thread
	for thread in results :

		#print "Reading %s" % thread

		# And over images in that thread
		for i in range(len(thread["name"])) :

			# If it's the right sheep, grab its results
			if thread["name"][i] == sheep :

				nonnormlac.append(thread["roylac"][i])
				entropy.append(thread["entropy"][i])
				entstd.append(thread["entstd"][i])
				directionality.append(thread["directionality"][i])
				scale.append(thread["scale"][i])
				id.append(thread["name"][i])
				imageid.append(thread["ID"][i])
				ratio.append(thread["ratio"][i])
				inflcount.append(thread["inflcount"][i])
				inflsize.append(thread["inflsize"][i])
				inflstd.append(thread["inflstd"][i])
				mindist.append(thread["mindist"][i])
				ifdist.append(thread["ifdist"][i])
				blur.append(thread["blur"][i])















# Write data to csv : Lacunarities, Entropy, Entropy std

nonnormlac = np.array(nonnormlac)
entropy = np.array(entropy)
entstd = np.array(entstd)

NLac = { "Sheep" : id, "Image" : imageid, "Lacunarity" : nonnormlac[:, 0].T }
Ent = { "Sheep" : id, "Image" : imageid, "Entropy" : entropy[:, 0].T }
Entstd = { "Sheep" : id, "Image" : imageid, "stdEntropy" : entstd[:, 0].T }



Dist = { "Sheep" : id, "Image" : imageid, "MinDist" : mindist, "IFDist" : ifdist }
Blur = { "Sheep" : id, "Image" : imageid, "Blur" : blur }



pd.DataFrame(NLac).to_csv("results/normalised_lacunarity.csv", columns = ["Sheep", "Image", "Lacunarity"])
pd.DataFrame(Ent).to_csv("results/entropy.csv", columns = ["Sheep", "Image", "Entropy"])
pd.DataFrame(Entstd).to_csv("results/entropy_std.csv", columns = ["Sheep", "Image", "stdEntropy"])
pd.DataFrame(Dist).to_csv("results/interfoci_dist.csv", columns = ["Sheep", "Image", "MinDist", "IFDist"])
pd.DataFrame(Blur).to_csv("results/blur.csv", columns = ["Sheep", "Image", "Blur"])





# Write data to csv : directionality, scale

Gabor = { "Sheep" : id, "Image" : imageid, "Scale" : np.array(scale).T, "Directionality" : np.array(directionality).T }

pd.DataFrame(Gabor).to_csv("results/gabor_filters.csv", columns = ["Sheep", "Image", "Scale", "Directionality"])







# Write data to csv : tissue to sinusoid ratio

Ratio = { "Sheep" : id, "Image" : imageid, "TSRatio" : np.array(ratio).T }
pd.DataFrame(Ratio).to_csv("results/tissue_sinusoid_ratio.csv", cols = ["Sheep", "Image", "TSRatio"])






# Write data to csv : inflammatory focus count and size

Foci = { "Sheep" : id, "Image" : imageid, "Count" : np.array(inflcount).T, "meanSize" : np.array(inflsize).T, "stdSize" : np.array(inflstd).T }
pd.DataFrame(Foci).to_csv("results/foci.csv", cols = ["Sheep", "Image", "Count", "meanSize", "stdSize"])







# Clean results directory
if not os.path.exists("results/threads") :
	os.makedirs("results/threads")

for old, new in zip(files, newfiles) :
	os.rename(old, new)






















"""
# Further : control vs Soay
control = {}
control["lacunarity"] = []
control["nonnormlac"] = []
control["entropy"] = []
control["entstd"] = []
control["directionality"] = []
control["scale"] = []

soay = copy.deepcopy(control)


for i, sheep in enumerate(id) :
	if sheep in [96, 97, 98, 99] : # if it's a control :
		print "CONTROL %d" % sheep
		control["lacunarity"].append(lacunarity[i, :])
		control["nonnormlac"].append(nonnormlac[i, :])
		control["entropy"].append(entropy[i, :])
		control["entstd"].append(entstd[i, :])
		control["directionality"].append(directionality[i])
		control["scale"].append(scale[i])

	else : # if it's a soay sheep
		soay["lacunarity"].append(lacunarity[i, :])
		soay["nonnormlac"].append(nonnormlac[i, :])
		soay["entropy"].append(entropy[i, :])
		soay["entstd"].append(entstd[i, :])
		soay["directionality"].append(directionality[i])
		soay["scale"].append(scale[i])
"""




















# Generate some plots - between soay and control

if not os.path.exists("results/plots") :
	os.makedirs("results/plots")


# Directionality
x = np.linspace(np.min(directionality), np.max(directionality), 200)[:, np.newaxis]
y = KernelDensity(bandwidth=0.03).fit(np.array(directionality)[:, np.newaxis])
#yc = KernelDensity(bandwidth=0.025).fit(np.array(control["directionality"])[:, np.newaxis])
#ys = KernelDensity(bandwidth=0.025).fit(np.array(soay["directionality"])[:, np.newaxis])


fig = plt.figure(figsize=(xdim*210/scalefactor, ydim*115/scalefactor))
plt.hist(directionality, 30, normed=True, color=colours[0], alpha=0.5)
plt.plot(x, np.exp(y.score_samples(x)), linewidth=3, c=colours[0])

#plt.hist(soay["directionality"], 20, normed=True, color=colours[2], alpha=0.5)
#plt.plot(x, np.exp(ys.score_samples(x)), linewidth=3, c=colours[2])

#plt.legend(["Control", "Soay"])
#plt.title("Directionality Coefficient Distribution\nKS Statistic = %.03f, p = %.2e" % ks_2samp(control["directionality"], soay["directionality"]))
plt.title("Directional Coefficients")
plt.xlabel("Directionality Coefficient")
plt.ylabel("Density")
plt.savefig("results/plots/directionality.pdf")
plt.close()

















# Scale
x = np.linspace(np.min(scale), np.max(scale), 200)[:, np.newaxis]
y = KernelDensity(bandwidth=1.2, kernel="gaussian").fit(np.array(scale)[:, np.newaxis])
#yc = KernelDensity(bandwidth=1.2, kernel="gaussian").fit(np.array(control["scale"])[:, np.newaxis])
#ys = KernelDensity(bandwidth=1.2, kernel="gaussian").fit(np.array(soay["scale"])[:, np.newaxis])


fig = plt.figure(figsize=(xdim*210/scalefactor, ydim*115/scalefactor))
plt.hist(scale, bins=np.unique(scale), normed=True, color=colours[0], alpha=0.5, align="left")
plt.plot(x, np.exp(y.score_samples(x)), linewidth=3, c=colours[0])
plt.xticks(np.unique(scale)[::2])

#plt.hist(soay["scale"], bins=bins, normed=True, color=colours[2], alpha=0.5)
#plt.plot(x, np.exp(ys.score_samples(x)), linewidth=3, c=colours[2])

#mwu, mwp = mannwhitneyu(control["scale"], soay["scale"])
#plt.legend(["Control", "Soay"])
plt.title("Characteristic Scales") # Distribution\nMann-Whitney U = %d, p = %.2e" % (mwu, 2*mwp))
plt.xlabel("Scale (pixels)")
plt.ylabel("Density")
plt.savefig("results/plots/scale.pdf")
plt.close()





















# Lacunarity and Entropy
# Plot at modal scales for both control and soay
#controlscale = mode(control["scale"])[0] * 4
#soayscale = mode(soay["scale"])[0] * 4

#controllac1 = np.array(control["lacunarity"])[:, np.where(lacscales == controlscale)[0]]
#controllac2 = np.array(control["lacunarity"])[:, np.where(lacscales == soayscale)[0]]

#soaylac1 = np.array(soay["lacunarity"])[:, np.where(lacscales == controlscale)[0]]
#soaylac2 = np.array(soay["lacunarity"])[:, np.where(lacscales == soayscale)[0]]


fig = plt.figure(figsize=(xdim*210/scalefactor, ydim*115/scalefactor))
#plt.subplot(121)

x = np.linspace(np.min(nonnormlac), np.max(nonnormlac), 200)[:, np.newaxis]
y = KernelDensity(bandwidth=0.0075).fit(np.array(nonnormlac[:, 0])[:, np.newaxis])
#yc = KernelDensity(bandwidth=0.015).fit(controllac1)
#ys = KernelDensity(bandwidth=0.015).fit(soaylac1)

plt.hist(nonnormlac, 30, normed=True, color=colours[0], alpha=0.5)
plt.plot(x, np.exp(y.score_samples(x)), linewidth=3, c=colours[0])

#plt.hist(soaylac1, 20, normed=True, color=colours[2], alpha=0.5)
#plt.plot(x, np.exp(ys.score_samples(x)), linewidth=3, c=colours[2])

#plt.legend(["Control", "Soay"])
plt.title("Normalised Lacunarity at Characteristic Scale") #\nKS Statistic = %.03f, p = %.2e" % ks_2samp(controllac1.squeeze(), soaylac1.squeeze()))
plt.xlabel("Normalised Lacunarity")
plt.ylabel("Density")

xl1 = np.where(y.score_samples(x) > -8)
xl2 = np.where(y.score_samples(x) > -8)
xrange = np.append(xl1, xl2)
plt.xlim(x[np.min(xrange)], x[np.max(xrange)])



"""
plt.subplot(122)

x = np.linspace(np.min(lacunarity), np.max(lacunarity), 200)[:, np.newaxis]
yc = KernelDensity(bandwidth=0.015).fit(controllac2)
ys = KernelDensity(bandwidth=0.015).fit(soaylac2)

plt.hist(controllac2, 20, normed=True, color=colours[0], alpha=0.5)
plt.plot(x, np.exp(yc.score_samples(x)), linewidth=3, c=colours[0])

plt.hist(soaylac2, 20, normed=True, color=colours[2], alpha=0.5)
plt.plot(x, np.exp(ys.score_samples(x)), linewidth=3, c=colours[2])

plt.legend(["Control", "Soay"])
plt.title("Normalised Lacunarity Distribution at Soay Scale\nKS Statistic = %.03f, p = %.2e" % ks_2samp(controllac2.squeeze(), soaylac2.squeeze()))
plt.xlabel("Normalised Lacunarity")
plt.ylabel("Density")

xl1 = np.where(yc.score_samples(x) > -8)
xl2 = np.where(ys.score_samples(x) > -8)
xrange = np.append(xl1, xl2)
plt.xlim(x[np.min(xrange)], x[np.max(xrange)])
"""
plt.tight_layout()

plt.savefig("results/plots/lacunarity.pdf")
plt.close()














"""

# Roy Lacunarity


# Lacunarity and Entropy
# Plot at modal scales for both control and soay
controlscale = mode(control["scale"])[0]
soayscale = mode(soay["scale"])[0]

controllac1 = np.array(control["nonnormlac"])[:, np.where(lacscales == controlscale * 2)[0]]
controllac2 = np.array(control["nonnormlac"])[:, np.where(lacscales == soayscale * 2)[0]]

soaylac1 = np.array(soay["nonnormlac"])[:, np.where(lacscales == controlscale * 2)[0]]
soaylac2 = np.array(soay["nonnormlac"])[:, np.where(lacscales == soayscale * 2)[0]]


fig = plt.figure(figsize=(xdim*210/scalefactor, ydim*115/scalefactor))
plt.subplot(121)

x = np.linspace(np.min(nonnormlac), np.max(nonnormlac), 200)[:, np.newaxis]
yc = KernelDensity(bandwidth=0.01).fit(controllac1)
ys = KernelDensity(bandwidth=0.01).fit(soaylac1)

plt.hist(controllac1, 20, normed=True, color=colours[0], alpha=0.5)
plt.plot(x, np.exp(yc.score_samples(x)), linewidth=3, c=colours[0])

plt.hist(soaylac1, 20, normed=True, color=colours[2], alpha=0.5)
plt.plot(x, np.exp(ys.score_samples(x)), linewidth=3, c=colours[2])

plt.legend(["Control", "Soay"])
plt.title("Normalised Lacunarity Distribution at Control Scale\nKS Statistic = %.03f, p = %.2e" % ks_2samp(controllac1.squeeze(), soaylac1.squeeze()))
plt.xlabel("Normalised Lacunarity")
plt.ylabel("Density")

xl1 = np.where(yc.score_samples(x) > -8)
xl2 = np.where(ys.score_samples(x) > -8)
xrange = np.append(xl1, xl2)
plt.xlim(x[np.min(xrange)], x[np.max(xrange)])




plt.subplot(122)

x = np.linspace(np.min(nonnormlac), np.max(nonnormlac), 200)[:, np.newaxis]
yc = KernelDensity(bandwidth=0.008).fit(controllac2)
ys = KernelDensity(bandwidth=0.008).fit(soaylac2)

plt.hist(controllac2, 20, normed=True, color=colours[0], alpha=0.5)
plt.plot(x, np.exp(yc.score_samples(x)), linewidth=3, c=colours[0])

plt.hist(soaylac2, 20, normed=True, color=colours[2], alpha=0.5)
plt.plot(x, np.exp(ys.score_samples(x)), linewidth=3, c=colours[2])

plt.legend(["Control", "Soay"])
plt.title("Normalised Lacunarity Distribution at Soay Scale\nKS Statistic = %.03f, p = %.2e" % ks_2samp(controllac2.squeeze(), soaylac2.squeeze()))
plt.xlabel("Normalised Lacunarity")
plt.ylabel("Density")

xl1 = np.where(yc.score_samples(x) > -8)
xl2 = np.where(ys.score_samples(x) > -8)
xrange = np.append(xl1, xl2)
plt.xlim(x[np.min(xrange)], x[np.max(xrange)])

plt.tight_layout()

plt.savefig("results/plots/roy_lac.pdf")
plt.close()

"""


















# Entropy


# Lacunarity and Entropy
# Plot at modal scales for both control and soay
#controlscale = mode(control["scale"])[0]
#soayscale = mode(soay["scale"])[0]

#controllac1 = np.array(control["entropy"])[:, np.where(lacscales == controlscale * 2)[0]]
#controllac2 = np.array(control["entropy"])[:, np.where(lacscales == soayscale * 2)[0]]

#soaylac1 = np.array(soay["entropy"])[:, np.where(lacscales == controlscale * 2)[0]]
#soaylac2 = np.array(soay["entropy"])[:, np.where(lacscales == soayscale * 2)[0]]


fig = plt.figure(figsize=(xdim*210/scalefactor, ydim*115/scalefactor))
#plt.subplot(121)

x = np.linspace(np.min(entropy), np.max(entropy), 200)[:, np.newaxis]
y = KernelDensity(bandwidth=0.005).fit(np.array(entropy[:, 0])[:, np.newaxis])
#ys = KernelDensity(bandwidth=0.005).fit(soaylac1)

plt.hist(entropy, 30, normed=True, color=colours[0], alpha=0.5)
plt.plot(x, np.exp(y.score_samples(x)), linewidth=3, c=colours[0])

#plt.hist(soaylac1, 20, normed=True, color=colours[2], alpha=0.5)
#plt.plot(x, np.exp(ys.score_samples(x)), linewidth=3, c=colours[2])

#plt.legend(["Control", "Soay"])
plt.title("Entropy at Characteristic Scale") #\nKS Statistic = %.03f, p = %.2e" % ks_2samp(controllac1.squeeze(), soaylac1.squeeze()))
plt.xlabel("Information Entropy")
plt.ylabel("Density")

xl1 = np.where(y.score_samples(x) > -8)
xl2 = np.where(y.score_samples(x) > -8)
xrange = np.append(xl1, xl2)
plt.xlim(x[np.min(xrange)], x[np.max(xrange)])



"""
plt.subplot(122)

x = np.linspace(np.min(entropy), np.max(entropy), 200)[:, np.newaxis]
yc = KernelDensity(bandwidth=0.004).fit(controllac2)
ys = KernelDensity(bandwidth=0.004).fit(soaylac2)

plt.hist(controllac2, 20, normed=True, color=colours[0], alpha=0.5)
plt.plot(x, np.exp(yc.score_samples(x)), linewidth=3, c=colours[0])

plt.hist(soaylac2, 20, normed=True, color=colours[2], alpha=0.5)
plt.plot(x, np.exp(ys.score_samples(x)), linewidth=3, c=colours[2])

plt.legend(["Control", "Soay"])
plt.title("Entropy Distribution at Soay Scale\nKS Statistic = %.03f, p = %.2e" % ks_2samp(controllac2.squeeze(), soaylac2.squeeze()))
plt.xlabel("Information Entropy")
plt.ylabel("Density")

xl1 = np.where(yc.score_samples(x) > -8)
xl2 = np.where(ys.score_samples(x) > -8)
xrange = np.append(xl1, xl2)
plt.xlim(x[np.min(xrange)], x[np.max(xrange)])
"""

plt.tight_layout()

plt.savefig("results/plots/entropy.pdf")
plt.close()


































# Inflammatory Foci
x = np.linspace(np.min(inflcount), np.max(inflcount), 200)[:, np.newaxis]
y = KernelDensity(bandwidth=1.2, kernel="gaussian").fit(np.array(inflcount)[:, np.newaxis])
#yc = KernelDensity(bandwidth=1.2, kernel="gaussian").fit(np.array(control["scale"])[:, np.newaxis])
#ys = KernelDensity(bandwidth=1.2, kernel="gaussian").fit(np.array(soay["scale"])[:, np.newaxis])


fig = plt.figure(figsize=(xdim*210/scalefactor, ydim*115/scalefactor))
plt.hist(inflcount, bins=np.unique(inflcount), normed=True, color=colours[0], alpha=0.5, align="left")
plt.plot(x, np.exp(y.score_samples(x)), linewidth=3, c=colours[0])

xl1 = np.where(y.score_samples(x) > -8)
xl2 = np.where(y.score_samples(x) > -8)
xrange = np.append(xl1, xl2)
plt.xlim(0, x[np.max(xrange)])

#plt.xticks([0, np.unique(scale)[::5])

#plt.hist(soay["scale"], bins=bins, normed=True, color=colours[2], alpha=0.5)
#plt.plot(x, np.exp(ys.score_samples(x)), linewidth=3, c=colours[2])

#mwu, mwp = mannwhitneyu(control["scale"], soay["scale"])
#plt.legend(["Control", "Soay"])
plt.title("Inflammatory Foci") # Distribution\nMann-Whitney U = %d, p = %.2e" % (mwu, 2*mwp))
plt.xlabel("Count")
plt.ylabel("Density")
plt.savefig("results/plots/foci.pdf")
plt.close()




















# Tissue to Sinusoid Ratio
x = np.linspace(np.min(ratio), np.max(ratio), 200)[:, np.newaxis]
y = KernelDensity(bandwidth=0.008, kernel="gaussian").fit(np.array(ratio)[:, np.newaxis])
#yc = KernelDensity(bandwidth=1.2, kernel="gaussian").fit(np.array(control["scale"])[:, np.newaxis])
#ys = KernelDensity(bandwidth=1.2, kernel="gaussian").fit(np.array(soay["scale"])[:, np.newaxis])


fig = plt.figure(figsize=(xdim*210/scalefactor, ydim*115/scalefactor))
plt.hist(ratio, 30, normed=True, color=colours[0], alpha=0.5)
plt.plot(x, np.exp(y.score_samples(x)), linewidth=3, c=colours[0])

xl1 = np.where(y.score_samples(x) > -8)
xl2 = np.where(y.score_samples(x) > -8)
xrange = np.append(xl1, xl2)
plt.xlim(x[np.min(xrange)], x[np.max(xrange)])

#plt.xticks([0, np.unique(scale)[::5])

#plt.hist(soay["scale"], bins=bins, normed=True, color=colours[2], alpha=0.5)
#plt.plot(x, np.exp(ys.score_samples(x)), linewidth=3, c=colours[2])

#mwu, mwp = mannwhitneyu(control["scale"], soay["scale"])
#plt.legend(["Control", "Soay"])
plt.title("Tissue to Sinusoid Ratio") # Distribution\nMann-Whitney U = %d, p = %.2e" % (mwu, 2*mwp))
plt.xlabel("Ratio")
plt.ylabel("Density")
plt.savefig("results/plots/ratio.pdf")
plt.close()











