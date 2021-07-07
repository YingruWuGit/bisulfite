# Clustering

This is an optional step. If the dataset contains sequences with different bisulfite accessible patches, then "clustering.py" is needed to divide the whole dataset into multiple groups.

In the "example" folder, the csv file "dsid_C_1461.csv" is a dataset needs to be divided into clusters. It's first column `id` is the id of each sequence and the last column `feq` is the frequency of each sequence. It's first row contains the location sites.

"clustering.py" needs two arguments:

Line 1 is the path and name of the csv file.

Line 2 is the cluster number.

It outputs a csv file for each cluster in the same folder.

# Bayesian Segmentation

"bayesian_segment.py" is used for detecting bisulfite accessible regions.

In the example folder, "dsid_C_1461_0.csv" is one of the clusters outputed by "clustering.py".

"bayesian_segment.py" needs only one argument, that is the path and name of dataset in line 1: ~/dsid_C_1461_0.csv

It outputs "dsid_C_1461_0_res.csv" in the same folder. The "dsid_C_1461_0_res.csv" will have two more rows.

Row 2 is the posterior probability of a change occured in each location.

Row 3 is the posterior probability that each location would respond to bisulfite.

It also outputs two graphs. One is a summary of dataset and posterior response rates, the other is the distribution of patch size.
