# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Robust Extraction of Quantitative Information from Histology Images

# <headingcell level=4>

# Quentin Caudron
# <br /><br />
# 
# Romain Garnier
# <br /><br />
# 
# *with Bryan Grenfell and Andrea Graham*

# <headingcell level=3>

# Outline

# <markdowncell>

# - Image processing
# - Extracted measures
# - Preliminary analysis
# - Future directions

# <markdowncell>

# 4. Age as random effect <---
# 
# ["interface_hepatitis", "confluent_necrosis", "portal_inflammation", "ln_ap_ri"]

# <codecell>

def normalise(df, skip = []) :
	for i in df.columns :
		if i not in skip :
			df[i] -= df[i].mean()
			df[i] /= df[i].std()
	return df






def rescale(df, skip = []) :
    for i in df.columns :
        if i not in skip :
            df[i] -= df[i].min()
            df[i] /= df[i].max()
    return df



# Remove a layer from a list
def delayer(m) :
	out = []
	for i in m :
		if isinstance(i, list) :
			for j in i :
				out.append(j)
		else :
			out.append(i)
	return out







# Remove all layers from a list
def flatten(m) :
	out = m[:]

	while out != delayer(out) :
		out = delayer(out)

	return out








# Generate all combinations of objects in a list
def combinatorial(l) :
	out = []

	for numel in range(len(l)) :
		for i in itertools.combinations(l, numel+1) :
			out.append(list(i))

	return out










def pcaplot(df) :

	# PCA
	pca = decomposition.PCA(whiten = True)
	pca.fit(df)
	p1 = pca.components_[0] / np.abs(pca.components_[0]).max() * np.sqrt(2)/2
	p2 = pca.components_[1] / np.abs(pca.components_[1]).max() * np.sqrt(2)/2

	# Normalise
	norms = np.max([np.sqrt((np.array(zip(p1, p2)[i])**2).sum()) for i in range(len(p1))])
	c = plt.Circle( (0, 0), radius = 1, alpha = 0.2)
	plt.axes(aspect = 1)
	plt.gca().add_artist(c)

	plt.scatter(p1 / norms, p2 / norms)
	plt.xlim([-1, 1])
	plt.ylim([-1, 1])

	for i, text in enumerate(df.columns) :
		plt.annotate(text, xy = [p1[i], p2[i]])

	plt.tight_layout()











def test_all_linear(df, y, x, return_significant = False, group = None) :

    # All possible combinations of independent variables
	independent = combinatorial(x)

	fits = {}
	pval = {}
	linmodels = {}
	qsum = {}
	aic = {}

	# For all dependent variables, one at a time
	for dependent in y :

		print "Fitting for %s." % dependent

		# For all combinations of independent variables
		for covariate in independent :

			# Standard mixed model
			if group is None :

				# Fit a linear model
				subset = delayer([covariate, dependent])
				df2 = df[delayer(subset)].dropna()
				df2["Intercept"] = np.ones(len(df2))
                
				ols = sm.GLS(endog = df2[dependent], exog = df2[delayer([covariate, "Intercept"])]).fit()

				# Save the results
				if (return_significant and ols.f_pvalue < 0.05) or (not return_significant) :
					linmodels.setdefault(dependent, []).append(ols)
					fits.setdefault(dependent, []).append(ols.rsquared)
					pval.setdefault(dependent, []).append(ols.f_pvalue)
					aic.setdefault(dependent, []).append(ols.aic)


			# Mixed effects model
			else :
				subset = delayer([covariate, dependent, group])
				df2 = df[delayer(subset)].dropna()

				# Fit a mixed effects model
				ols = MixedLM(endog = df2[dependent], exog = df2[covariate], groups = df2[group]).fit()

				# Calculate AIC
				linmodels.setdefault(dependent, []).append(ols)
				fits.setdefault(dependent, []).append(2 * (ols.k_fe + 1) - 2 * ols.llf)
				pval.setdefault(dependent, []).append(ols.pvalues)

	if group is not None :
		for i in y :
			f = np.array(fits[i])
			models = np.array(linmodels[i])
			idx = np.where(f - f.min() <= 2)[0]
			bestmodelDoF = [j.k_fe for j in np.array(linmodels[i])[idx]]
			bestmodels = [idx[j] for j in np.where(bestmodelDoF == np.min(bestmodelDoF))[0]]
			qsum[i] = models[idx[np.where(f[bestmodels] == np.min(f[bestmodels]))]]


		return linmodels, fits, pval, qsum

	return linmodels, fits, pval, aic

	
		















def summary(models) :

	# Generate list of everything
	r2 = np.array([m.r2 for dependent in models.keys() for m in models[dependent]])
	p = np.array([m.f_stat["p-value"] for dependent in models.keys() for m in models[dependent]])
	mod = np.array([m for dependent in models.keys() for m in models[dependent]])
	dependent = np.array([dependent for dependent in models.keys() for m in models[dependent]])

	# Sort by R2
	idx = np.argsort(r2)[::-1]

	# Output string
	s = "%d significant regressions.\n\n" % len(r2)
	s += "Ten most correlated :\n\n"

	# Print a summary of the top ten correlations
	for i in idx[:10] :
		s += ("%s ~ %s\n" % (dependent[i], " + ".join(mod[i].x.columns[:-1])))
		s += ("R^2 = %f\tp = %f\n\n" % (r2[i], p[i]))

	print s
    
    
    
    
    
    
    
def rstr(y, x) :
    formatstr = "%s ~ " % y
    for i in x[:-1] :
        formatstr += str(i)
        formatstr += " + "
    formatstr += str(x[-1])
    return formatstr








# <codecell>

import numpy as np
from sklearn.neighbors import KernelDensity
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn
import pandas as pd
import itertools
from sklearn import linear_model, ensemble, decomposition, cross_validation, preprocessing
from statsmodels.regression.mixed_linear_model import MixedLM
import statsmodels.api as sm
from statsmodels.regression.linear_model import OLSResults
from statsmodels.tools.tools import add_constant


%matplotlib inline
rcParams["figure.figsize"] = (14, 8)


# RAW DATA

raw_physical = pd.read_csv("../data/physical.csv")
raw_histo = pd.read_csv("../data/tawfik.csv")
ent = pd.read_csv("../4x/results/entropy.csv").drop(["Unnamed: 0"], 1)
foci = pd.read_csv("../4x/results/foci.csv").drop(["Unnamed: 0"], 1)
lac = pd.read_csv("../4x/results/normalised_lacunarity.csv").drop(["Unnamed: 0"], 1)
gabor = pd.read_csv("../4x/results/gabor_filters.csv").drop(["Unnamed: 0"], 1)
ts = pd.read_csv("../4x/results/tissue_sinusoid_ratio.csv").drop(["Unnamed: 0"], 1)

raw_image = pd.merge(lac, ent,
	on=["Sheep", "Image"]).merge(foci, 
	on=["Sheep", "Image"]).merge(gabor,
	on=["Sheep", "Image"]).merge(ts, 
    on=["Sheep", "Image"])
raw_image.rename(columns = {	"meanSize" : "FociSize", 
								"TSRatio" : "TissueToSinusoid",
								"Count" : "FociCount" }, inplace=True)



# CLEAN DATA

physcols = ["Weight", "Sex", "AgeAtDeath", "Foreleg", "Hindleg"]
imagecols = ["Entropy", "Lacunarity", "Inflammation", "Scale", "Directionality", "FociCount", "FociSize", "TissueToSinusoid"]
histcols = ["Lobular_collapse", "Interface_hepatitis", "Confluent_necrosis", "Ln_ap_ri", "Portal_inflammation", "BD_hyperplasia", "Fibrosis", "TawfikTotal", "Mean_hep_size", "Min_hep_size", "Max_hep_size"]





# IMAGE

# Set FociSize to zero if FociCount is zero
# Drop stdSize
image = raw_image
image = image.drop("stdSize", 1)
image.FociSize[raw_image.FociCount == 0] = 0



# HISTO

histo = raw_histo
histo = histo.drop(["Vessels", "Vacuol", "Pigment", "Std_hep_size"], 1)



# PHYSICAL

physical = raw_physical
physical = physical.drop(["CurrTag", "DeathDate", "Category"], 1)
physical




# COMPLETE DATASET

raw_data = pd.merge(pd.merge(image, histo, on="Sheep", how="outer"), physical, on="Sheep", how="outer")
raw_data.to_csv("../data/tentative_complete.csv")




# AVERAGED BY SHEEP
data = raw_data
data["Inflammation"] = data.FociCount * data.FociSize

sheep = rescale(data.groupby("Sheep").mean())
age = rescale(data.groupby("AgeAtDeath").mean())







# REGRESSIONS : fixed effects, grouped by sheep

df = sheep[["Portal_inflammation", "FociSize"]].dropna()
df["Intercept"] = np.ones(len(df))
portal_inflammation = sm.GLS(endog = df.Portal_inflammation, exog = df[["FociSize", "Intercept"]]).fit().summary()
#portal_inflammation = portal_inflammation.summary()
del portal_inflammation.tables[2]



df = sheep[["BD_hyperplasia", "Scale", "Directionality", "FociSize"]].dropna()
df["Intercept"] = np.ones(len(df))
hyperplasia = sm.GLS(endog = df.BD_hyperplasia, exog = df[["FociSize", "Scale", "Directionality", "Intercept"]]).fit().summary()
#hyperplasia.summary()
del hyperplasia.tables[2]






# REGRESSIONS : fixed effects, grouped by age

df = age[["Max_hep_size", "Entropy", "Directionality"]].dropna()
df["Intercept"] = np.ones(len(df))
maxhepsize = sm.GLS(endog = df.Max_hep_size, exog = df[["Entropy", "Directionality", "Intercept"]]).fit().summary()
del maxhepsize.tables[2]




df = age[["Lobular_collapse", "FociSize"]].dropna()
df["Intercept"] = np.ones(len(df))
lobular_collapse = sm.GLS(endog = df.Lobular_collapse, exog = df[["FociSize", "Intercept"]]).fit().summary()
del lobular_collapse.tables[2]


df = age[["Interface_hepatitis", "Lacunarity"]].dropna()
df["Intercept"] = np.ones(len(df))
interface_hepatitis = sm.GLS(endog = df.Interface_hepatitis, exog = df[["Lacunarity", "Intercept"]]).fit().summary()
del interface_hepatitis.tables[2]


df = age[["Fibrosis", "Inflammation"]].dropna()
df["Intercept"] = np.ones(len(df))
fibrosis = sm.GLS(endog = df.Fibrosis, exog = df[["Inflammation", "Intercept"]]).fit().summary()
del fibrosis.tables[2]




# PCA

s = sheep.dropna(subset=delayer([imagecols, histcols]))
pca = decomposition.PCA(n_components=1)
pcax = pca.fit_transform(s[imagecols])
pcay = pca.fit_transform(s[histcols])
pca = sm.GLS(endog = pcay[:, 0][:, np.newaxis], exog = add_constant(pcax)).fit().summary()
del pca.tables[2]





# REGRESSIONS : mixed effects, intercept on age at death

df = age[["Fibrosis", "Inflammation"]].dropna()
df["Intercept"] = np.ones(len(df))
fibrosis = sm.GLS(endog = df.Fibrosis, exog = df[["Inflammation", "Intercept"]]).fit().summary()
del fibrosis.tables[2]

# <codecell>

a = portal_inflammation.summary()
del a.tables[2]
a

# <headingcell level=2>

# Image Processing

# <markdowncell>

# <img src="figures/sheep.jpg"></img>

# <markdowncell>

# <img src="figures/processed.jpg"></img>

# <headingcell level=3>

# Extraction

# <markdowncell>

# - Automagical
# - Reasonably quick

# <headingcell level=3>

# Robust

# <markdowncell>

# - Invariant to staining, slicing, field-related variation
# - Capture intersample variation

# <markdowncell>

# ![image](figures/robust3.jpg)

# <markdowncell>

# ![image](figures/robust4.jpg)

# <markdowncell>

# ![image](figures/robust1.jpg)

# <markdowncell>

# ![image](figures/robust2.jpg)

# <headingcell level=2>

# Structural and Textural Measures

# <markdowncell>

# - characteristic **scale** of sinusoid widths
# - **directional** amplitude of preferred sinusoid alignment
# - **tissue to sinusoid** ratio
# - **count** of inflammatory foci per image
# - **mean size** of inflammatory foci per image
# - information **entropy** of sinusoid distribution
# - **lacunarity** ( clustering ) of sinusoids

# <markdowncell>

# <img src="figures/gif.gif"></img>

# <markdowncell>

# ![image](figures/intra.png)

# <markdowncell>

# ![image](figures/inter2.png)

# <headingcell level=2>

# Exploratory Analysis

# <headingcell level=3>

# by individual

# <codecell>

portal_inflammation

# <markdowncell>

# ![image](figures/portal_inflammation.png)

# <codecell>

hyperplasia

# <markdowncell>

# ![image](figures/hyperplasia.png)

# <codecell>

pca

# <markdowncell>

# ![image](figures/pca.png)

# <headingcell level=2>

# Exploratory Analysis

# <headingcell level=3>

# by age class

# <codecell>

fibrosis

# <markdowncell>

# ![image](figures/fibrosis.png)

# <codecell>

lobular_collapse

# <markdowncell>

# ![image](figures/lobular_collapse.png)

# <codecell>

interface_hepatitis

# <markdowncell>

# ![image](figures/interface_hepatitis.png)

# <headingcell level=2>

# Exploratory analysis

# <headingcell level=3>

# with a random effect on age at death

# <markdowncell>

# | Dependent variable                       | Models<br />AIC < 2 + AIC<sub>min</sub> | Primary explanatory variables                           |
# |------------------------------------------|:----------------------------------:|---------------------------------------------------------------------|
# | Ishak score                              |                  7                 | entropy, tissue-to-sinusoid, focus count, focus size                |
# | Lobular collapse                         |                  5                 | entropy, lacunarity, tissue-to-sinusoid, focus count                |
# | Confluent necrosis                       |                  1                 | entropy                                                             |
# | Interface hepatitis                      |                  2                 | entropy, tissue-to-sinusoid                                         |
# | Portal inflammation                      |                  4                 | entropy, focus size, lacunarity, focus count, scale, directionality |
# | Fibrosis                                 |                  2                 | entropy, lacunarity, tissue-to-sinusoid                             |
# | Biliary hyperplasia                      |                  1                 | focus size                                                          |
# | Necrosis, apoptosis, random inflammation |    <font color="white">This_is_bla</font>2<font color="white">This_is_bla</font>           | entropy, lacunarity                                                 |

# <markdowncell>

# - entropy consistently explains histological measures when controlled for age
# - also important : tissue to sinusoid ratio, focus count and size, lacunarity

# <markdowncell>

# - biological / historical reasoning for this potential cohort effect
# - interpretation of these models
# - quality of fit

# <headingcell level=2>

# Conclusions

# <markdowncell>

# - our **semi-educated guess** measures may capture relevant information
# - underlying **structure** in the data needs thought
# - still no **map** from image or histological measures to condition of individual

# <headingcell level=2>

# Future directions

# <headingcell level=3>

# Further exploration of the dataset

# <markdowncell>

# - 145 sheep ( 89 females )
# - 11 age classes
# - potential redundancy in various measures

# <markdowncell>

# - 4460 entries across 27 variables
# - 3330 with full image and histological information
# - 1196 for which **complete** information is available

# <headingcell level=3>

# More data

# <markdowncell>

# - nutritional information
# - immunity data

# <headingcell level=3>

# Narrow-field images

# <markdowncell>

# - 12536 images
# - spatial distribution of nuclei

# <markdowncell>

# ![image](figures/10.jpg)

# <markdowncell>

# ![image](figures/Processed2.jpg)

# <markdowncell>

# ![image](figures/Segmented.jpg)

# <markdowncell>

# <img src="figures/10x.png" width=100%></src>

