# Clustering

This is an optional step. If the dataset contains sequences with different bisulfite accessible patches, then "clustering.py" is needed to divide the whole dataset into multiple groups.

In the "example" folder, the csv file "dsid_C_1461.csv" is a dataset needs to be divided into clusters. It's first column `id` is the id of each sequence and the last column `feq` is the frequency of each sequence. It's first row contains the location sites.

"clustering.py" needs two arguments:

Line 1 is the path and name of the csv file.

Line 2 is the cluster number.

It outputs a csv file for each cluster with the same format, in the same folder.

# Bayesian Segmentation

"bayesian_segment.py" is used for detecting bisulfite accessible regions (BAR).

In the "example" folder, "dsid_C_1461_0.csv" is one of the clusters outputed by "clustering.py".

"bayesian_segment.py" needs only one argument:

Line 1 is the path and name of the csv file.

It outputs "dsid_C_1461_0_res.csv" in the same folder. The "dsid_C_1461_0_res.csv" will have two additional rows:

Row 2 `prob of change` is the posterior probability of a change occured for each site.

Row 3 `prob of response` is the posterior response probability to bisulfite for each site.

It also outputs two graphs. One is the summary of dataset and posterior response rates, the other is the distribution of patch size.

![alt text](https://github.com/YingruWuGit/bisulfite/blob/main/example/Figure_1.png)

In the above fig, the horizontal axis is th bp locations and the vertical axis is the response probability to bisulfite.

The blue dots and lines are the means and standard deviations form the data, the red dots are the posterior estimations form our Bayesian segmentation model.

![alt text](https://github.com/YingruWuGit/bisulfite/blob/main/example/Figure_2.png=250x250)

In the above fig, the blue bars are the patch size distribution.
