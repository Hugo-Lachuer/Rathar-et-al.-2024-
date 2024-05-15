# Rathar-et-al.-2024-
Datasets and codes for spatial analysis

1/ DatasetPillar
This folder contains the data from cells seeded on pillars. The folder contains 3 types of files (below X is the ID of the cell):
-“P-image X.csv”: This file contains the x and y coordinates of all the actin structures.
-“P-image X_grid.csv”: This file contains the properties of the pillars grid. Images have been rotated to have pillar orientations parallel to the x and y axes (i.e. crystal directions [10] and [01] parallel to x and y axes). Four pillars are measured for each field of view. The “Feret” column gives the measured diameter of the pillar while the “X” and “Y” columns give the coordinates. The pillar grid can be extended from these 4 pillars.
-“P-Image X_mask.tif” : Image file giving the cell mask.
In addition, the folder contains a “ImageDimensions.txt” which gives, for each image, the size in pixels and µm to obtain the pixel size.

2/ DatasetFlat
This folder contains the data from cells seeded on a flat surface. It contains the same types of files as “DatasetPillar” except for the pillar coordinates files, which are not defined for a flat surface.

3/ Spatial analysis of actin structures Dataset1 (pillar)
This code imports “datasetPillar” and computes the Ripley’s K functions and its multivariate form (i.e. K-cross function) as well as the associated CSR simulations.

4/ Spatial analysis of actin structures Dataset2 (flat)
This code imports “datasetFlat” and computes the Ripley’s K functions as well as the associated CSR simulations. Moreover, the code imports results from both pillar and flat datasets, plots them together and compares Ripley’s K functions using a permutation test.

5/ Functions
This .R file contains functions called by “Spatial analysis of actin structures Dataset1 (pillar)” and “Spatial analysis of actin structures Dataset2 (flat)”.
