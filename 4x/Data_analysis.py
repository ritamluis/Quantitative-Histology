# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn
import pandas as pd
import itertools
import os
from sklearn import linear_model, ensemble, decomposition, cross_validation, preprocessing
from statsmodels.regression.mixed_linear_model import MixedLM
import statsmodels
import statsmodels.api as sm
from statsmodels.regression.linear_model import OLSResults
from statsmodels.tools.tools import add_constant
from sklearn.neighbors import KernelDensity

from mpl_toolkits.mplot3d import Axes3D, proj3d

print statsmodels.__version__

%matplotlib inline
rcParams["figure.figsize"] = (14, 8)

rcParams["text.usetex"] = False
rcParams["xtick.labelsize"] = 12
rcParams["ytick.labelsize"] = 12
rcParams["font.size"] = 14
rcParams["axes.titlesize"] = 16
#rcParams["text.usetex"] = False
rcParams["font.family"] = "Serif"
rcParams["figure.dpi"] = 600

# <codecell>

# RAW DATA

raw_physical = pd.read_csv("../data/physical.csv")
raw_histo = pd.read_csv("../data/tawfik.csv")
ent = pd.read_csv("results/entropy.csv").drop(["Unnamed: 0"], 1)
foci = pd.read_csv("results/foci.csv").drop(["Unnamed: 0"], 1)
lac = pd.read_csv("results/normalised_lacunarity.csv").drop(["Unnamed: 0"], 1)
gabor = pd.read_csv("results/gabor_filters.csv").drop(["Unnamed: 0"], 1)
ts = pd.read_csv("results/tissue_sinusoid_ratio.csv").drop(["Unnamed: 0"], 1)
blur = pd.read_csv("results/blur.csv").drop(["Unnamed: 0"], 1)
distances = pd.read_csv("results/interfoci_dist.csv").drop(["Unnamed: 0"], 1)

raw_image = pd.merge(lac, ent,
	on=["Sheep", "Image"]).merge(foci, 
	on=["Sheep", "Image"]).merge(gabor,
	on=["Sheep", "Image"]).merge(ts, 
	on=["Sheep", "Image"]).merge(blur, 
	on=["Sheep", "Image"]).merge(distances, 
    on=["Sheep", "Image"])
raw_image.rename(columns = {	"meanSize" : "FociSize", 
								"TSRatio" : "TissueToSinusoid",
								"Count" : "FociCount" }, inplace=True)

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
def combinatorial(l, short = np.inf) :
	out = []

	for numel in range(len(l)) :
		for i in itertools.combinations(l, numel+1) :
			if len(list(i)) < short :
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









def big_ass_matrix(df, y, x, group = None, short = True) :

    independent = combinatorial(x, short)
    
    models = {}
    p = {}
    aic = {}
    r2 = {}
    best = {}
    dfs = {}
    bestdf = {}
    
    for dependent in y :
        
        print "Regressing for %s" % dependent
        
        for covariate in independent :
            
            if group is None :
                
                subset = delayer([covariate, dependent])
                df2 = df[subset].dropna()
                df2["Intercept"] = np.ones(len(df2))
                dfs.setdefault(dependent, []).append(df2)
                
                ols = sm.GLS(endog=df2[dependent], exog=df2[delayer([covariate, "Intercept"])]).fit()
                
                models.setdefault(dependent, []).append(ols)
                p.setdefault(dependent, []).append(ols.pvalues[:-1].values)
                aic.setdefault(dependent, []).append(ols.aic)
                r2.setdefault(dependent, []).append(ols.rsquared)
            
            else :
                
                subset = delayer([covariate, dependent, group])
                df2 = df[subset].dropna()
                dfs.setdefault(dependent, []).append(df2)
                
                ols = MixedLM.from_formula(rstr(y=dependent, x=covariate), data=df2, groups=df2[group]).fit()
                
                models.setdefault(dependent, []).append(ols)
                aic.setdefault(dependent, []).append(2 * (ols.k_fe + 1) - 2 * ols.llf)
                p.setdefault(dependent, []).append(ols.pvalues[1:-1].values)
                r2.setdefault(dependent, []).append(mmR2(df2, ols))

    
       
        bestAIC = np.min(aic[dependent])
        
        for i, val in enumerate(models[dependent]) :
            
            if aic[dependent][i] < 2 + bestAIC :
                
                if np.sum(p[dependent][i] > 0.05) == 0 :
                    
                    if group is None :
                        
                        best.setdefault(dependent, []).append(val)
                        bestdf.setdefault(dependent, []).append(dfs[dependent][i])
                        
                    else :
                        
                        if val.random_effects.abs().mean()[0] > 0.01 :
                            
                            best.setdefault(dependent, []).append(val)
                            bestdf.setdefault(dependent, []).append(dfs[dependent][i])
       
            
        if best.has_key(dependent) :
            for i, model in enumerate(best[dependent]) :
                
                if not os.path.exists("regressions/%s" % dependent) :                  
                    os.mkdir("regressions/%s" % dependent)
                    
                if not os.path.exists("../talk/figures/regressions/%s" % dependent) :                  
                    os.mkdir("../talk/figures/regressions/%s" % dependent)                

                if group is None :

                    dfx = bestdf[dependent][i]
                    plt.scatter(model.fittedvalues.values, dfx[model.model.endog_names].values, c=seaborn.color_palette("deep", 8)[0])
                    plt.plot(dfx[model.model.endog_names].values, dfx[model.model.endog_names].values, c=seaborn.color_palette("deep", 8)[2])
                    plt.ylabel(model.model.endog_names)
                    yl = model.model.exog_names[:]
                    yl.remove("Intercept")
                    plt.xlabel("Estimate using " + ", ".join(yl))
                    plt.title(rstr(dependent, model.model.exog_names).replace(" + Intercept", ""))
                    #plt.title(r"$R^2$ = %.02f" % model.rsquared)
                    st = ("$R^2$ = %.03f\n\n"% model.rsquared)
                    for coefnum, coef in enumerate(yl) :
                        st += ("%s" % coef)
                        st += (" : %.03f\n" % model.params[coef])
                        st += ("$p$ = %.01e\n\n" % model.pvalues[coefnum])
                    #plt.suptitle(st)
                    plt.text(0.01, .99, st, va="top", ha="left")
                    plt.xlim([-0.05, 1.05])
                    plt.ylim([-0.05, 1.05])
                    plt.savefig("regressions/%s/lm-%d.pdf" % (dependent, i))
                    plt.savefig("../talk/figures/regressions/%s/lm-%d.png" % (dependent, i), dpi=300, jpeg_quality=90)
                    plt.close()

                else :
                    dfx = bestdf[dependent][i]
                    y, yhat = mmPredict(model.model.data.frame, model)
                    plt.scatter(yhat, y, c=seaborn.color_palette("deep", 8)[0])
                    plt.plot(y, y, c=seaborn.color_palette("deep", 8)[2])
                    plt.ylabel(model.model.endog_names)
                    yl = model.model.exog_names[:]
                    yl.remove("Intercept")
                    plt.xlabel("Estimate using " + ", ".join(yl))
                    plt.title(rstr(dependent, model.model.exog_names).replace("Intercept + ", ""))
                    
                    #plt.title(r"$R^2$ = %.02f" % mmR2(dfx, model))
                    st = ("$R^2$ = %.03f\n\n" % mmR2(dfx, model))
                    for coefnum, coef in enumerate(yl) :
                        st += coef
                        st += " : %.03f\n" % model.fe_params[1+coefnum]
                        st += "$p$ = %.01e\n\n" % model.pvalues[coef]
                    st += ("Avg. abs. RE coef. : %.03f" % model.random_effects.abs().mean())
                    plt.text(0.01, .99, st, va="top", ha="left")
                    
                    plt.xlim([-0.05, 1.05])
                    plt.ylim([-0.05, 1.05])
                    plt.savefig("regressions/%s/mm_%d.pdf" % (dependent, i))
                    plt.savefig("../talk/figures/regressions/%s/mm_%d.png" % (dependent, i), dpi=300, jpeg_quality=90)
                    plt.close()
        
    return best, (models, p, r2, aic)
    
    
    
    














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









def mmR2(df, ols) :

    y = df[ols.model.endog_names]
    sigma_a = ols.random_effects.values.var()

    yhat = np.zeros_like(y)
    for i, coef in enumerate(ols.fe_params) :
        yhat += ols.model.data.exog[:, i] * coef

    sigma_f = yhat.var()

    idxx = []
    
    grouplabel = list(df.columns)
    grouplabel.remove(ols.model.endog_names)
    remove_these = ols.model.exog_names[:]
    remove_these.remove("Intercept")
    for i in remove_these :
        grouplabel.remove(i)

    for j, u in enumerate(ols.model.group_labels) :
        idx = np.where(df[grouplabel] == u)[0]
        yhat[idx] += ols.random_effects.values[j]

    sigma_e = np.var(yhat - y)

    return (sigma_f + sigma_a) / (sigma_f + sigma_a + sigma_e)





def mmPredict(df, ols) :

    y = df[ols.model.endog_names]
    sigma_a = ols.random_effects.values.var()

    yhat = np.zeros_like(y)
    for i, coef in enumerate(ols.fe_params) :
        yhat += ols.model.data.exog[:, i] * coef

    sigma_f = yhat.var()

    idxx = []
    
    grouplabel = list(df.columns)
    grouplabel.remove(ols.model.endog_names)
    remove_these = ols.model.exog_names[:]
    remove_these.remove("Intercept")
    for i in remove_these :
        grouplabel.remove(i)

    for j, u in enumerate(ols.model.group_labels) :
        idx = np.where(df[grouplabel] == u)[0]
        yhat[idx] += ols.random_effects.values[j]
        
    return y, yhat




#best, stuff = big_ass_matrix(df=sheep, y=["Weight"], x=imagecols, group="AgeAtDeath", short=True)

# <codecell>

raw_image.columns

# <codecell>

# CLEAN DATA

physcols = ["Weight", "Sex", "AgeAtDeath", "Foreleg", "Hindleg", "E", "CES", "RES"]
imagecols = ["Entropy", "Lacunarity", "Inflammation", "Scale", "Directionality", "FociCount", "FociSize", "TissueToSinusoid", "Blur", "MinDist", "IFDist"]
histcols = ["LobularCollapse", "InterfaceHepatitis", "ConfluentNecrosis", "LnApRi", "PortalInflammation", "BDHyperplasia", "Fibrosis", "TawfikTotal", "MeanHepSize", "MinHepSize", "MaxHepSize"]





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
exposure = pd.read_csv("../data/exposure.csv")
exposure.rename(columns={"BirthYear" : "BirthYear", "AvgOfLambWS" : "E"}, inplace=True)
exposure["AgeAtDeath"] = 2011 - exposure.BirthYear
exposure.E = np.round(1. - exposure.E, 6)
exposure["CES"] = np.flipud(np.flipud(exposure.E).cumsum())
exposure["RES"] = np.flipud(np.flipud(exposure.E - exposure.CES / (exposure.AgeAtDeath + 1)).cumsum())
physical = pd.merge(physical, exposure, on="AgeAtDeath", how="inner")




# COMPLETE DATASET

raw_data = pd.merge(pd.merge(image, histo, on="Sheep", how="outer"), physical, on="Sheep", how="outer")
raw_data.to_csv("../data/tentative_complete.csv")




# AVERAGED BY SHEEP
data = raw_data
data["Inflammation"] = data.FociCount * data.FociSize

sheep = np.round(rescale(data.groupby("Sheep").mean()), 9)



sheep.rename(columns = {"Lobular_collapse" : "LobularCollapse", 
                        "Interface_hepatitis" : "InterfaceHepatitis",
                        "Confluent_necrosis" : "ConfluentNecrosis",
                        "Ln_ap_ri" : "LnApRi",
                        "Portal_inflammation" : "PortalInflammation",
                        "BD_hyperplasia" : "BDHyperplasia",
                        "Mean_hep_size" : "MeanHepSize",
                        "Min_hep_size" : "MinHepSize",
                        "Max_hep_size" : "MaxHepSize"}, inplace=True)



data.rename(columns = {"Lobular_collapse" : "LobularCollapse", 
                        "Interface_hepatitis" : "InterfaceHepatitis",
                        "Confluent_necrosis" : "ConfluentNecrosis",
                        "Ln_ap_ri" : "LnApRi",
                        "Portal_inflammation" : "PortalInflammation",
                        "BD_hyperplasia" : "BDHyperplasia",
                        "Mean_hep_size" : "MeanHepSize",
                        "Min_hep_size" : "MinHepSize",
                        "Max_hep_size" : "MaxHepSize"}, inplace=True)

#age = rescale(data.groupby("Sheep").mean().dropna(subset=delayer([imagecols, histcols, "AgeAtDeath"])).groupby("AgeAtDeath").mean())

# <codecell>

pcols = physcols[:]
pcols.remove("AgeAtDeath")

print " -------- 1"
best_lm_phys, stuff_lm_phys = big_ass_matrix(df=sheep, y=physcols, x=imagecols, group=None, short=5)

print " -------- 2"
best_lm_hist, stuff_lm_hist = big_ass_matrix(df=sheep, y=histcols, x=imagecols, group=None, short=5)

print " -------- 3"
best_mm_phys, stuff_mm_phys = big_ass_matrix(df=sheep, y=pcols, x=imagecols, group="AgeAtDeath", short=5)

print " -------- 4"
best_mm_hist, stuff_mm_hist = big_ass_matrix(df=sheep, y=histcols, x=imagecols, group="AgeAtDeath", short=5)

# <codecell>

y = "BDHyperplasia"
x = ["Inflammation", "Scale", "Directionality"]

dfx = sheep[delayer([x, y, "AgeAtDeath"])].dropna()
model = MixedLM.from_formula(rstr(y, x), data=dfx, groups="AgeAtDeath").fit()
#model = sm.GLS(endog=dfx.Portal_inflammation, exog=dfx[["FociSize", "AgeAtDeath"]]).fit()

dfx = sheep[["BDHyperplasia", "Inflammation", "AgeAtDeath"]].dropna()
model2 = MixedLM.from_formula(rstr(y, ["Inflammation"]), data=dfx, groups="AgeAtDeath").fit()

dfx = sheep[["BDHyperplasia", "FociSize", "AgeAtDeath"]].dropna()
model3 = MixedLM.from_formula(rstr(y, ["FociSize"]), data=dfx, groups="AgeAtDeath").fit()

# <codecell>

ss = "E"
s = np.array([sheep[sheep.AgeAtDeath == model.random_effects.index.values[i]][ss].iloc[0] for i in range(len(model.random_effects.index.values))])
s -= s.min()
s /= s.max()

idx = np.round(s * 2000).astype(int)
C = seaborn.color_palette("coolwarm", 2001)
CC = []
for i in idx :
    CC.append(C[i])


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.scatter(model.random_effects, model2.random_effects, zs=model3.random_effects, s=450, c=CC, alpha=1)#, c=np.array(C)[idx.astype(int)], alpha=1)
plt.title("Biliary Hyperplasia : Coefficients of Birth Year Random Effect\nFirst component of PCA explains 98.3%% of variance.\nWarmer colours - higher %s" % ss, y=0.9)

for i, sheep_age in enumerate(model.random_effects.index) :
    sheep_e = sheep[sheep.AgeAtDeath == sheep_age].E.iloc[0]
    x2, y2, _ = proj3d.proj_transform(model.random_effects.values[i], model2.random_effects.values[i], model3.random_effects.values[i], ax.get_proj())
    plt.annotate(int(np.round(2011 - sheep_age * 11)), xy=(x2, y2), xytext=(x2-0.0105, y2-0.002))

ax.set_axis_bgcolor("white")
plt.savefig("../talk/figures/regressions/BDHyperplasia/mm_coefs_color_%s.png" % ss, dpi=600, jpeg_quality=100)

# <codecell>

plt.scatter(exposure.BirthYear, exposure.CES)

# <codecell>

plt.scatter(model.random_effects.index.values, idx)

# <codecell>


#plt.plot([sheep[np.round(sheep.AgeAtDeath * 11) == i].E.iloc[0] for i in range(11)])
#plt.plot(exposure.AgeAtDeath, exposure.E)
C2 = []
for i in (np.unique(sheep[["AgeAtDeath", "BD_hyperplasia"]].dropna().AgeAtDeath) * 11).astype(int) :
    C2.append(sheep[np.round(sheep.AgeAtDeath * 11) == i].E.iloc[0])
#sheep[["AgeAtDeath", "E"]].dropna()

# <codecell>

plt.plot(exposure.BirthYear, exposure.E)

# <codecell>

exposure

# <codecell>

plt.text()

# <codecell>

df = sheep[["Interface_hepatitis", "AgeAtDeath", "Entropy"]].dropna()
ols = MixedLM.from_formula(rstr("Interface_hepatitis", ["Entropy"]), data=df, groups=df.AgeAtDeath).fit()
ols.summary()

# <codecell>

print mmR2(df, ols)
print ols.random_effects.abs().mean()

# <codecell>

print ols.summary()

print ols.bse
print ols.bse_re
print ols.bse_fe

# <codecell>

histo.columns

# <codecell>

y1, y2 = mmPredict(df, ols)
plt.scatter(y1, y2)
plt.plot(y1,y1)

# <codecell>


# <codecell>

ols = pd.stats.api.ols(y = avesc.Portal_inflammation, x = avesc.Inflammation)
print ols.summary

# <codecell>

plt.scatter(avesc.FociSize, avesc.Portal_inflammation)
plt.plot(avesc.FociSize, avesc.FociSize * ols.beta[0] + ols.beta[1])
plt.title("p = %.02e" % ols.p_value[0])

# <codecell>

ols2 = pd.stats.api.ols(y = avesc.BD_hyperplasia, x = avesc[["FociSize", "Directionality", "Scale"]])
print ols2.summary

# <codecell>

X = ["FociSize", "Directionality", "Scale"]
x = (ols2.beta[-1] + 
    ols2.beta[0] * avesc[X[0]] +
    ols2.beta[1] * avesc[X[1]] + 
    ols2.beta[2] * avesc[X[2]])
y = avesc.BD_hyperplasia

plt.scatter(x, y)
#plt.plot([0, 1], [- 0.0504, 1- 0.0504] )
#plt.plot(y, x)

#plt.scatter(avesc.FociSize, avesc.Portal_inflammation)
#plt.plot(avesc.FociSize, avesc.FociSize * ols.beta[0] + ols.beta[1])
#plt.title("p = %.02e" % ols.p_value[0])

# <codecell>

ols3 = pd.stats.api.ols(y = avesc.Portal_inflammation, x = avesc.Inflammation)
plt.scatter(avesc.Inflammation, avesc.Portal_inflammation)
plt.plot(avesc.Inflammation, avesc.Inflammation * ols3.beta[0] + ols3.beta[1])
print ols3.summary

# <codecell>

ols4 = pd.stats.api.ols(y = avesc.QTotal, x = avesc.TawfikTotal)
plt.scatter(avesc.TawfikTotal, avesc.QTotal)
plt.plot(avesc.TawfikTotal, avesc.TawfikTotal * ols4.beta[0] + ols4.beta[1])
print ols4.summary

# <codecell>

df2 = sheep[["TawfikTotal", "Entropy", "AgeAtDeath"]]
df2.dropna(inplace=True)
ols = MixedLM.from_formula(rstr("TawfikTotal", ["Entropy"]), data=df2, groups="AgeAtDeath").fit()

ols.pvalues[1:-1].values

# <codecell>

df = sheep[["Inflammation", "Mean_hep_size", "FociSize"]].dropna()
df["Intercept"] = np.ones(len(df))
ols = sm.GLS(endog=df.Mean_hep_size, exog=df[["Inflammation", "FociSize"]]).fit()
ols.rsquared

# <codecell>

mmR2(df, ols)






print "YO"

# <codecell>

print "YO"

# <codecell>

y = pca.fit_transform(s[histcols])

# <codecell>

mm, mmfits, mmpvals, mmqsum = test_all_linear(sheep, ["TawfikTotal"], imagecols, group="AgeAtDeath")

# <codecell>

df = sheep[delayer(["TawfikTotal", "Inflammation", "ResidualES"])].dropna()
df["Intercept"] = np.ones(len(df))
tt = MixedLM(endog = df.TawfikTotal, exog = df[["Inflammation"]], groups=df.ResidualES).fit()
tt.summary()
#del fibrosis.tables[2]


# <codecell>

models, fits, pvals, blah = test_all_linear(sheep, ["Ln_ap_ri"], imagecols, group="AgeAtDeath")

# <codecell>

np.where(np.array(fits["Ln_ap_ri"]) < (2 + np.min(fits["Ln_ap_ri"])))

# <codecell>

idx = 0
print models["Ln_ap_ri"][idx].summary()
print fits["Ln_ap_ri"][idx]

# <codecell>

a.model.endog_names

plt.scatter(y, yhat, c=seaborn.color_palette("deep", 8)[0])
plt.plot(y, y, c=seaborn.color_palette("deep", 8)[2])
plt.xlabel(a.model.endog_names)
yl = a.model.exog_names
#yl.remove("Intercept")
plt.ylabel(", ".join(yl))
plt.title("R2 = %.02f")
st = "%s : %.03f, p = %.03e.\n" * len(yl)
stl = []
for i in range(len(yl)) :
    stl.append(yl[i])
    stl.append(p[dependent][i])
plt.suptitle(st % tuple(delayer([yl[i], p[dependent][i] for i in range(len(yl))

# <codecell>

a.model.exog_names

# <codecell>


# <codecell>

raw_data[raw_data.TissueToSinusoid == raw_data.TissueToSinusoid.max()]

# <codecell>

df = sheep[["MeanHepSize", "Directionality"]].dropna()
df["Intercept"] = np.ones(len(df))
ols = sm.GLS(endog=df.MeanHepSize, exog=df[["Directionality", "Intercept"]]).fit()
ols.summary()

# <codecell>


