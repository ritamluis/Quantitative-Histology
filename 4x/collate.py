# Collate data from numerous batch.py runs 
# into one environment and into .csv files

import copy
import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn
from scipy.stats import mode, ks_2samp, mannwhitneyu
from sklearn.neighbors import KernelDensity





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

lacunarity = []
nonnormlac = []
entropy = []
entstd = []
directionality = []
imageid = []
id = []
scale = []


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

				lacunarity.append(thread["lacunarity"][i])
				nonnormlac.append(thread["roylac"][i])
				entropy.append(thread["entropy"][i])
				entstd.append(thread["entstd"][i])
				directionality.append(thread["directionality"][i])
				scale.append(thread["scale"][i])
				id.append(thread["name"][i])
				imageid.append(thread["ID"][i])






# Write data to csv : Lacunarities, Entropy, Entropy std
gaborscales = 1. / np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]) # 2 to 20 pixels
lacscales = np.unique(np.append(np.round(np.exp(np.linspace(np.log(4), np.log(400), 15))), 2/gaborscales)).astype(int)

lacunarity = np.array(lacunarity)
nonnormlac = np.array(nonnormlac)
entropy = np.array(entropy)
entstd = np.array(entstd)

order = np.insert(lacscales.astype(str), 0, ["Sheep", "Image"])

Lac = { "Sheep" : id, "Image" : imageid }
NLac = { "Sheep" : id, "Image" : imageid }
Ent = { "Sheep" : id, "Image" : imageid }
Entstd = { "Sheep" : id, "Image" : imageid }

for i, s in enumerate(lacscales) :
	Lac["%d" % s] = lacunarity[:, i]
	NLac["%d" % s] = nonnormlac[:, i]
	Ent["%d" % s] = entropy[:, i]
	Entstd["%d" % s] = entstd[:, i]


pd.DataFrame(Lac).to_csv("results/normalised_lacunarity.csv", cols = order)
pd.DataFrame(NLac).to_csv("results/raw_lacunarity.csv", cols = order)
pd.DataFrame(Ent).to_csv("results/entropy.csv", cols = order)



# Write data to csv : directionality, scale

Gabor = { "Sheep" : id, "Image" : imageid, "Scale" : np.array(scale).T, "Directionality" : np.array(directionality).T }

pd.DataFrame(Gabor).to_csv("results/gabor_filters.csv", cols = ["Sheep", "Image", "Scale", "Directionality"])







# Clean results directory
if not os.path.exists("results/threads") :
	os.makedirs("results/threads")

for old, new in zip(files, newfiles) :
	os.rename(old, new)




















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





















# Generate some plots - between soay and control

if not os.path.exists("results/plots") :
	os.makedirs("results/plots")


# Directionality
x = np.linspace(np.min(directionality), np.max(directionality), 200)[:, np.newaxis]
yc = KernelDensity(bandwidth=0.025).fit(np.array(control["directionality"])[:, np.newaxis])
ys = KernelDensity(bandwidth=0.025).fit(np.array(soay["directionality"])[:, np.newaxis])

plt.hist(control["directionality"], 20, normed=True, color=colours[0], alpha=0.5)
plt.plot(x, np.exp(yc.score_samples(x)), linewidth=3, c=colours[0])

plt.hist(soay["directionality"], 20, normed=True, color=colours[2], alpha=0.5)
plt.plot(x, np.exp(ys.score_samples(x)), linewidth=3, c=colours[2])

plt.legend(["Control", "Soay"])
plt.title("Directionality Coefficient Distribution\nKS Statistic = %.03f, p = %.2e" % ks_2samp(control["directionality"], soay["directionality"]))
plt.xlabel("Directionality Coefficient")
plt.ylabel("Density")
plt.savefig("results/plots/directionality.pdf")
plt.close()

















# Scale
x = np.linspace(np.min(scale), np.max(scale), 200)[:, np.newaxis]
yc = KernelDensity(bandwidth=1.2, kernel="gaussian").fit(np.array(control["scale"])[:, np.newaxis])
ys = KernelDensity(bandwidth=1.2, kernel="gaussian").fit(np.array(soay["scale"])[:, np.newaxis])

bins = [4.5, 5.5, 6.5, 7.5, 9, 11, 13.5, 18]

plt.hist(control["scale"], bins=bins, normed=True, color=colours[0], alpha=0.5)
plt.plot(x, np.exp(yc.score_samples(x)), linewidth=3, c=colours[0])

plt.hist(soay["scale"], bins=bins, normed=True, color=colours[2], alpha=0.5)
plt.plot(x, np.exp(ys.score_samples(x)), linewidth=3, c=colours[2])

mwu, mwp = mannwhitneyu(control["scale"], soay["scale"])
plt.legend(["Control", "Soay"])
plt.title("Characteristic Scale Distribution\nMann-Whitney U = %d, p = %.2e" % (mwu, 2*mwp))
plt.xlabel("Scale (pixels)")
plt.ylabel("Density")
plt.savefig("results/plots/scale.pdf")
plt.close()





















# Lacunarity and Entropy
# Plot at modal scales for both control and soay
controlscale = mode(control["scale"])[0]
soayscale = mode(soay["scale"])[0]

controllac1 = np.array(control["lacunarity"])[:, np.where(lacscales == controlscale * 2)[0]]
controllac2 = np.array(control["lacunarity"])[:, np.where(lacscales == soayscale * 2)[0]]

soaylac1 = np.array(soay["lacunarity"])[:, np.where(lacscales == controlscale * 2)[0]]
soaylac2 = np.array(soay["lacunarity"])[:, np.where(lacscales == soayscale * 2)[0]]


plt.subplot(121)

x = np.linspace(np.min(lacunarity), np.max(lacunarity), 200)[:, np.newaxis]
yc = KernelDensity(bandwidth=0.015).fit(controllac1)
ys = KernelDensity(bandwidth=0.015).fit(soaylac1)

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

plt.tight_layout()

plt.savefig("results/plots/lacunarity.pdf")
plt.close()















# Roy Lacunarity


# Lacunarity and Entropy
# Plot at modal scales for both control and soay
controlscale = mode(control["scale"])[0]
soayscale = mode(soay["scale"])[0]

controllac1 = np.array(control["nonnormlac"])[:, np.where(lacscales == controlscale * 2)[0]]
controllac2 = np.array(control["nonnormlac"])[:, np.where(lacscales == soayscale * 2)[0]]

soaylac1 = np.array(soay["nonnormlac"])[:, np.where(lacscales == controlscale * 2)[0]]
soaylac2 = np.array(soay["nonnormlac"])[:, np.where(lacscales == soayscale * 2)[0]]


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




















# Entropy


# Lacunarity and Entropy
# Plot at modal scales for both control and soay
controlscale = mode(control["scale"])[0]
soayscale = mode(soay["scale"])[0]

controllac1 = np.array(control["entropy"])[:, np.where(lacscales == controlscale * 2)[0]]
controllac2 = np.array(control["entropy"])[:, np.where(lacscales == soayscale * 2)[0]]

soaylac1 = np.array(soay["entropy"])[:, np.where(lacscales == controlscale * 2)[0]]
soaylac2 = np.array(soay["entropy"])[:, np.where(lacscales == soayscale * 2)[0]]


plt.subplot(121)

x = np.linspace(np.min(entropy), np.max(entropy), 200)[:, np.newaxis]
yc = KernelDensity(bandwidth=0.005).fit(controllac1)
ys = KernelDensity(bandwidth=0.005).fit(soaylac1)

plt.hist(controllac1, 20, normed=True, color=colours[0], alpha=0.5)
plt.plot(x, np.exp(yc.score_samples(x)), linewidth=3, c=colours[0])

plt.hist(soaylac1, 20, normed=True, color=colours[2], alpha=0.5)
plt.plot(x, np.exp(ys.score_samples(x)), linewidth=3, c=colours[2])

plt.legend(["Control", "Soay"])
plt.title("Entropy Distribution at Control Scale\nKS Statistic = %.03f, p = %.2e" % ks_2samp(controllac1.squeeze(), soaylac1.squeeze()))
plt.xlabel("Information Entropy")
plt.ylabel("Density")

xl1 = np.where(yc.score_samples(x) > -8)
xl2 = np.where(ys.score_samples(x) > -8)
xrange = np.append(xl1, xl2)
plt.xlim(x[np.min(xrange)], x[np.max(xrange)])




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

plt.tight_layout()

plt.savefig("results/plots/entropy.pdf")
plt.close()








