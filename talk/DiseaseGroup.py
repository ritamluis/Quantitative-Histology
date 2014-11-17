# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn
from matplotlib import rcParams
rcParams["figure.figsize"] = (14, 8)
rcParams["xtick.labelsize"] = 12
rcParams["ytick.labelsize"] = 12
rcParams["font.size"] = 14
rcParams["axes.titlesize"] = 16
#rcParams["text.usetex"] = False
rcParams["font.family"] = "Serif"
rcParams["figure.dpi"] = 600


a = pd.read_csv("../data/villagebay_population.csv")
b = pd.read_csv("../data/exposure.csv")

fig, (ax, ax2) = plt.subplots(2, 1, sharex=True)

#ax = plt.subplot(211)
ax.plot(a.Year, a.VillageBay, c=seaborn.color_palette("deep", 8)[0], lw=3)
ax.scatter(a.Year, a.VillageBay, c=seaborn.color_palette("deep", 8)[0], s=50)
ax.set_title("Village Bay Population")
ax.set_ylim([180, 700])


#ax2 = plt.subplot(212, sharex=ax)

ax2.plot(b.BirthYear, b.AvgOfLambWS, c=seaborn.color_palette("deep", 8)[2], lw=3)
ax2.scatter(b.BirthYear, b.AvgOfLambWS, c=seaborn.color_palette("deep", 8)[2], s=50)
ax2.set_title("Lamb Winter Survival")
ax2.set_xlim([1984.5, 2013.5])
ax2.set_ylim([0, 0.8])

plt.savefig("figures/population2.jpg", dpi=300, jpeg_quality=100)

# <headingcell level=1>

# Robust Extraction of Quantitative Information from Histology Images

# <headingcell level=4>

# Quentin Caudron

# <headingcell level=2>

# The Soay Sheep

# <markdowncell>

# <img src="figures/graphics/soay.jpg" />

# <markdowncell>

# <img src="figures/graphics/population2.jpg" />

# <markdowncell>

# <img src="figures/graphics/lit1.jpg" />

# <markdowncell>

# <img src="figures/graphics/lit2.jpg" />

# <markdowncell>

# <img src="figures/graphics/lit4.jpg" />

# <headingcell level=2>

# Outline

# <markdowncell>

# - Methods and data collection
# - Image processing
# - Extracted measures
# - Preliminary analysis
# - Future directions

# <headingcell level=2>

# Data

# <markdowncell>

# **In the field, winter of 2011 - 2012 :**
#     
# - Daily study area monitoring for deaths
# - 143 liver samples collected within a day of death

# <markdowncell>

# **In the lab :**
# 
# - Sectioning after paraffin treatment
# - H&E staining of about 1000 slides

# <markdowncell>

# **Analysis :**
# 
# - Pathology standard : semi-quantitative scoring
# - Image processing

# <headingcell level=3>

# The Field &copy;

# <markdowncell>

# Sweat-and-blood-collected in cold, cold Scotland.

# <markdowncell>

# Eight physical measurements :
# - Age at death
# - Weight
# - Sex
# - Limb length
# - Environmental "stress"

# <headingcell level=3>

# Clinical Pathology

# <markdowncell>

# Operator-driven visual analysis of 98 slides under microscopy.

# <markdowncell>

# Eleven discrete and continuous measures :
# 
# - Inflammation
# - Necrosis
# - Apoptosis
# - Hyperplasia
# - Fibrosis
# - Hepatitis

# <headingcell level=3>

# Image Processing

# <markdowncell>

# Automated analysis of 4430 images of slides representing 143 sheep.

# <markdowncell>

# Seven structural and textural measures with varying levels of biological interpretation :
# 
# - Inflammation
# - Hyperplasia / tissue density
# - Best-guess proxies for "generic degeneration"

# <headingcell level=2>

# Image Processing

# <markdowncell>

# <img src="figures/graphics/sheep.jpg"></img>

# <markdowncell>

# <img src="figures/graphics/processed.jpg"></img>

# <headingcell level=3>

# The Challenge

# <markdowncell>

# **Information extraction must be**
# - automagical - no operator input
# - reasonably quick - restricted computing time
# - robust - invariant to slicing, staining, field-related variation 
# - unbiased - same algorithms for everyone

# <markdowncell>

# ![image](figures/graphics/robust3.jpg)

# <markdowncell>

# ![image](figures/graphics/robust4.jpg)

# <markdowncell>

# ![image](figures/graphics/robust1.jpg)

# <markdowncell>

# ![image](figures/graphics/robust2.jpg)

# <markdowncell>

# <img src="figures/graphics/gif.gif"></img>

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

# ![image](figures/graphics/intra.png)

# <markdowncell>

# ![image](figures/graphics/inter2.png)

# <headingcell level=2>

# Exploratory Analysis

# <headingcell level=3>

# by individual

# <markdowncell>

# <img src="figures/regressions/BDHyperplasia/lm-0.png" />

# <markdowncell>

# <img src="figures/regressions/PortalInflammation/lm-0.png" />

# <markdowncell>

# <img src="figures/regressions/PortalInflammation/lm-1.png" />

# <headingcell level=2>

# Exploratory Analysis

# <headingcell level=3>

# controlled for age / cohort

# <markdowncell>

# <img src="figures/regressions/PortalInflammation/mm_0.png" />

# <markdowncell>

# <img src="figures/regressions/BDHyperplasia/mm_0.png" />

# <markdowncell>

# <img src="figures/regressions/BDHyperplasia/mm_1.png" />

# <markdowncell>

# <img src="figures/regressions/TawfikTotal/mm_0.png" />

# <markdowncell>

# <img src="figures/regressions/Fibrosis/mm_0.png" />

# <markdowncell>

# <img src="figures/regressions/PortalInflammation/mm_0.png" />

# <markdowncell>

# <img src="figures/regressions/Hindleg/mm_0.png" />

# <markdowncell>

# <img src="figures/regressions/Weight/mm_0.png" />

# <headingcell level=2>

# Further analysis

# <headingcell level=3>

# Age or cohort effect ?

# <markdowncell>

# <img src="figures/regressions/BDHyperplasia/mm_coefs_color_E.png" />

# <markdowncell>

# <img src="figures/regressions/BDHyperplasia/mm_coefs_color_CES.png" />

# <markdowncell>

# <img src="figures/regressions/BDHyperplasia/mm_coefs_color_RES.png" />

# <headingcell level=2>

# Conclusions

# <markdowncell>

# - our image measures capture **relevant** and **useful** information
# - a number of correlations can be **explained** biologically
# - underlying **structure** in the data needs thought
# - still no **map** from image or histological measures to condition of individual

# <headingcell level=2>

# Future directions

# <headingcell level=3>

# Further exploration of the dataset

# <markdowncell>

# - 145 sheep ( 89 females )
# - 12 age classes
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

# ![image](figures/graphics/10.jpg)

# <markdowncell>

# ![image](figures/graphics/Processed2.jpg)

# <markdowncell>

# ![image](figures/graphics/Segmented.jpg)

# <markdowncell>

# <img src="figures/graphics/10x.png" width=100%></src>

# <headingcell level=2>

# With thanks to

# <markdowncell>

# Romain Garnier
# 
# Andrea Graham
# 
# Tawfik Aboellail (CSU)
# 
# Bryan Grenfell

