# Clustering

This is an optional step. If the dataset contains sequences with different bisulfite accessible patches, then "clustering.py" is used to divide the whole dataset into multiple groups.

In the example folder, "dsid_C_1461.csv" is a dataset needs clustering. It's first column is the id of each sequence and the last column is the frequency of each sequence.

"clustering.py" needs two arguments.

In line 1, input the path and name of the file: ~/dsid_C_1461.csv

In line 2, input the cluster number.

It outputs multiple csv files for the clusters in the same folder.

# Bayesian Segmentation

"bayesian_segment.py" is used for detecting bisulfite accessible regions.

In the example folder, "dsid_C_1461_0.csv" is one of the clusters outputed by "clustering.py".

"bayesian_segment.py" needs only one argument, that is the path and name of dataset in line 1: ~/dsid_C_1461_0.csv

It outputs "dsid_C_1461_0_res.csv" in the same folder. The "dsid_C_1461_0_res.csv" will have two more rows.

Row 2 is the posterior probability of a change occured in each location.

Row 3 is the posterior probability that each location would respond to bisulfite.

It also outputs two graphs. One is a summary of dataset and posterior response rates, the other is the distribution of patch size.
