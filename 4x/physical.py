import numpy as np
import matplotlib.pyplot as plt
import seaborn
import pandas as pd
import itertools
from sklearn import linear_model, ensemble, decomposition, cross_validation, preprocessing
from statsmodels.regression.mixed_linear_model import MixedLM









def normalise(df, skip = []) :
	for i in df.columns :
		if i not in skip :
			df[i] -= df[i].mean()
			df[i] /= df[i].std()
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











def test_all_linear(df, y, x, return_significant = True, group = None) :

	# All possible combinations of independent variables
	independent = combinatorial(x)

	fits = {}
	pval = {}
	linmodels = {}
	qsum = {}

	# For all dependent variables, one at a time
	for dependent in y :

		print "Fitting for %s." % dependent

		# For all combinations of independent variables
		for covariate in independent :

			# Standard mixed model
			if group is None :

				# Fit a linear model
				ols = pd.stats.api.ols(y = df[dependent], x = df[covariate])

				# Save the results
				if (return_significant and ols.f_stat["p-value"] < 0.05) or (not return_significant) :
					linmodels.setdefault(dependent, []).append(ols)
					fits.setdefault(dependent, []).append(ols.r2)
					pval.setdefault(dependent, []).append(ols.f_stat["p-value"])


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

	return linmodels, fits, pval

	
		















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
	#return s

























# Read physical data
physical = pd.read_csv("../physical.csv").drop(["CurrTag", "DeathDate", "Category"], 1)



# Read improc data
ent = pd.read_csv("results/entropy.csv").drop(["Unnamed: 0"], 1)
foci = pd.read_csv("results/foci.csv").drop(["Unnamed: 0"], 1)
lac = pd.read_csv("results/normalised_lacunarity.csv").drop(["Unnamed: 0"], 1)
gabor = pd.read_csv("results/gabor_filters.csv").drop(["Unnamed: 0"], 1)
ts = pd.read_csv("results/tissue_sinusoid_ratio.csv").drop(["Unnamed: 0"], 1)



# Merge dataframes
sheep = np.unique(ent.Sheep)

Ent = [np.mean(ent.loc[ent.Sheep == i]).Entropy for i in sheep]
EntVar = [np.var(ent.loc[ent.Sheep == i]).Entropy for i in sheep]

Focicount = [np.mean(foci.loc[foci.Sheep == i]).Count for i in sheep]
FocicountVar = [np.var(foci.loc[foci.Sheep == i]).Count for i in sheep]
Focisize = [np.mean(foci.loc[foci.Sheep == i]).meanSize for i in sheep]
FocisizeVar = [np.var(foci.loc[foci.Sheep == i]).meanSize for i in sheep]

Lac = [np.mean(lac.loc[lac.Sheep == i]).Lacunarity for i in sheep]
LacVar = [np.var(lac.loc[lac.Sheep == i]).Lacunarity for i in sheep]

Scale = [np.mean(gabor.loc[gabor.Sheep == i]).Scale for i in sheep]
ScaleVar = [np.var(gabor.loc[gabor.Sheep == i]).Scale for i in sheep]
Dir = [np.mean(gabor.loc[gabor.Sheep == i]).Directionality for i in sheep]
DirVar = [np.var(gabor.loc[gabor.Sheep == i]).Directionality for i in sheep]

TS = [np.mean(ts.loc[ts.Sheep == i]).TSRatio for i in sheep]
TSVar = [np.var(ts.loc[ts.Sheep == i]).TSRatio for i in sheep]

improc = pd.DataFrame({	"Sheep" : sheep,
						"Entropy" : Ent,
						"EntropyVar" : EntVar,
						"FociCount" : Focicount,
						"FociCountVar" : FocicountVar,
						"FociSize" : Focisize,
						"FociSizeVar" : FocisizeVar,
						"Lacunarity" : Lac,
						"LacunarityVar" : LacVar,
						"Scale" : Scale,
						"ScaleVar" : ScaleVar,
						"Directionality" : Dir,
						"DirectionalityVar" : DirVar,
						"TissueToSinusoid" : TS,
						"TissueToSinusoidVar" : TSVar
						})




physcols = ["Weight", "Sex", "AgeAtDeath", "Foreleg", "Hindleg"]
imagecols = ["Entropy", "Lacunarity", "Scale", "Directionality", "FociCount", "FociSize", "TissueToSinusoid"]




# Merges :

# Sheep-centred dataframe
rawsheepdata = pd.merge(physical, improc, on="Sheep", how="outer")
sheepdata = normalise(rawsheepdata, skip = "Sheep")

# Image-centred dataframe
rawimagedata = pd.merge(lac, ent,
	on=["Sheep", "Image"]).merge(foci.drop("stdSize", 1), 
	on=["Sheep", "Image"]).merge(gabor,
	on=["Sheep", "Image"]).merge(ts, 
	on=["Sheep", "Image"]).merge(normalise(physical, skip = "Sheep"),
	on="Sheep")
rawimagedata.rename(columns = {	"meanSize" : "FociSize", 
								"TSRatio" : "TissueToSinusoid",
								"Count" : "FociCount" }, inplace=True)
imagedata = normalise(rawimagedata, skip = physcols + ["Sheep"])






# COLUMN ON NUMBER OF IMAGES






















"""
# CV

m = [	[[linear_model.ElasticNet(max_iter=1000000, alpha = i, l1_ratio = j) 
			for i in np.linspace(0.1, 1.0, 100)]
			for j in np.linspace(0., 1.0, 100)],
		[[[[linear_model.BayesianRidge(n_iter = 10000, alpha_1 = a1, alpha_2 = a2, lambda_1 = l1, lambda_2 = l2) 
			for a1 in np.logspace(1e-8, 1e-4, 20)] 
			for a2 in np.logspace(1e-8, 1e-4, 20)]
			for l1 in np.logspace(1e-8, 1e-4, 20)]
			for l2 in np.logspace(1e-8, 1e-4, 20)],
		[ensemble.RandomForestRegressor(n_estimators = i) for i in np.unique(np.logspace(1, 4, 100).astype(int))],
		[[[[ensemble.GradientBoostingRegressor(n_estimators = e, loss=l, learning_rate = r, max_depth = m) 
			for l in ["ls", "lad", "huber"]] 
			for e in np.unique(np.logspace(1, 4, 100).astype(int))]
			for r in np.logspace(-4, 0, 200)]
			for m in np.arange(1, 20)],
		#ensemble.AdaBoostRegressor(n_estimators = 100),
		#ensemble.AdaBoostRegressor(base_estimator = ensemble.GradientBoostingRegressor(n_estimators = 100), n_estimators = 100),
		#ensemble.AdaBoostRegressor(base_estimator = ensemble.RandomForestRegressor(n_estimators = 50), n_estimators = 100),
		[ensemble.BaggingRegressor(n_estimators = e) for e in np.unique(np.logspace(1, 4, 100).astype(int))]
	]

methods = flatten(m)





# For each physical measurement, perform fits using all above methods
scores = {}



for col in physcols :
	for q, method in enumerate(methods[10000:10100]) :
		if not scores.has_key(col) :
			scores[col] = []
		scores[col].append(
			cross_validation.cross_val_score(
				method, d[imagecols], d[col], cv=cross_validation.ShuffleSplit(len(d), 100, test_size=0.1), 
			n_jobs = -1, scoring = "r2").mean())
		print "Done %d" % q

"""



