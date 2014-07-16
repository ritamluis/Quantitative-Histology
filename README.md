Quantitative Histology
======================

Context
-------

Tissue samples can contain significant information on the state of health of an individual, and their 
interpretation can aid in the elucidation of mechanisms of degeneration of health. Commonly, a small sample of 
images is captured from histological preparations, and these images are assessed by eye. This can be 
prohibitatively slow at the population level, and human error can cause bias in the results. In addition, the 
procurement of quantitative information allows the use of rigorous statistical methods which would otherwise 
be unavailable to the analysis.

This repository contains in-progress code for the structural analysis of histopathological slides of soay 
sheep liver. Analyses occur on two levels of zoom in images. 

Widefield : 2.0 x 3.0 mm
------------------

Representative of the general state of the tissue. This is wide enough to see a number of inflammatory foci 
and a signficant amount of hepatic tissue.

- Characteristic scale of sinusoids ( by Gabor filtering )
- Directionality coefficient ( by Gabor filtering )
- Information entropy at characteristic scale
- Lacunarity at characteristic scale
- Count of inflammatory foci
- Approximate average size of inflammatory foci
- Ratio of tissue to sinusoids ( a measure of tissue density )

Narrowfield : 0.8 x 1.2 mm
--------------------------

Images are targetted at hepatic portal triads. 

- Background and inflamed nuclear density
- Count of nuclei ( watershed segmentation for overlapping nuclei )
- Nuclear radial distribution function
- Dimensions and areas of inflamed portal triads
