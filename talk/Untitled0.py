# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Robust Extraction of Quantitative Information from Histology Images

# <markdowncell>

# **Quentin Caudron**

# <markdowncell>

# <img src="figures/graphics/soay.jpg" />

# <markdowncell>

# <img src="figures/population3.png" width=1200px/>

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

# <br />**In the lab :**
# 
# - Sectioning after paraffin treatment
# - H&E staining of about 1000 slides

# <markdowncell>

# <br />**Analysis :**
# 
# - Pathology standard : semi-quantitative scoring
# - Image processing

# <headingcell level=3>

# The Field

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

# <img src="figures/graphics/gif.gif"></img>

# <markdowncell>

# <img src="figures/graphics/intra3.png" width=100%/>

# <markdowncell>

# <img src="figures/graphics/inter3.png"/>

# <headingcell level=2>

# Exploratory Analysis

# <markdowncell>

# by individual

# <markdowncell>

# <img src="figures/regressions/BDHyperplasia/lm-0.png" />

# <markdowncell>

# <img src="figures/regressions/PortalInflammation/lm-0.png" />

# <markdowncell>

# <img src="figures/regressions/PortalInflammation/lm-1.png" />

# <headingcell level=2>

# Exploratory Analysis

# <markdowncell>

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

# <img src="figures/regressions/Hindleg/mm_0.png" />

# <markdowncell>

# <img src="figures/regressions/Weight/mm_0.png" />

# <headingcell level=2>

# Further analysis

# <markdowncell>

# Age or cohort effect ?

# <markdowncell>

# <img src="figures/regressions/BDHyperplasia/new/mm_coefs_color_E.png" width=85%/>

# <markdowncell>

# <img src="figures/regressions/BDHyperplasia/new/mm_coefs_color_CES.png" width=85%/>

# <markdowncell>

# <img src="figures/regressions/BDHyperplasia/new/mm_coefs_color_RES.png" width=85%/>

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

# - 4460 entries across 31 variables
# - 3596 with full image and histological information
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

