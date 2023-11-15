This is a subset of routines from IRI tool helps for conducting kmeans analysis.
The bulk of this code is a wrapper around call to kmeans(), in order to compute
the "classifiability index" from Michealangelo along with the clustering output.

Main call: [CI,K]=kmeans_ci(X,standardize_flag,percent_variance_to_retain,clusters,repeat)
X[n,p] is time x grids
CI is classifiability index
K[n,1] is cluster assignment
This is just a loop around call to matlab function kmeans()

The utility kmeans_ar_ci_test() gives the CIs for random red noise, kind of like a confidence interval
Choose the first cluster number that is outside of the confidence interval.

Alternate routine kmeans_ci2/kmeans_ar_ci_test2 -- does an internal loop on a range of cluster sizes to give CIs
Also takes the average of the CI for a cluster, not the minimum



L Agel changes:

USE:  kmeans_ci_LA and kmeans_ar_ci_test_LA!!!  

This does not use the loop, but better because I save each CI and red-noise
CI matrix for each cluster size -- you can quit the program if taking too 
long and restart -- will not destroy earlier work.

Other important differences:
1. Changes kmeans_ci to calculate CI more like Michelangeli
2. Does all preprocessing outside of kmeans_ci (standardizing, removing daily 
mean, weighting for latitude, etc), so it is only done once.
3. If the data is serially correlated, calls ebisuzaki instead of kmeans_ar_ci_test 
4. Makes sure to seed the random generator properly.