An implementation of Bayesian Switching Dynamical Systems (BSDS)

If you found it useful in your work, please cite: 
Taghia J, Cai W, Ryali S, et al. Uncovering hidden brain state dynamics that regulate performance and decision-making during cognition. Nature Communications. 2018;9:2505. doi:10.1038/s41467-018-04723-6.


*** Usage ***
The main function is < BayesianSwitchingDynamicalSystems.m >. The default setting works best for fewer number of ROIs. By default, max_ldim is set to one less than the number of ROIs, that is max_ldim=n_roi-1. In applications that you have large number of ROIs, you may want to use the advanced initialization. In particular, you would need to bound max_ldim, for example to account for some variations in the data (say ~90%).

Let's say you have trained a group-level BSDS model. Of course, from the group model you know the estimated group covariance matrices. Now, if you wish to also know the subject-level covariances, you can use < compute_subject_level_stats.m > which takes as one of its input arguments, the group-level trained model, and then outputs the estimated covariances for each subject.

If you have questions about the usage and related issues, please write us an email: jalil.taghia@it.uu.se
wdcai@stanford.edu

Happy computing!
JAlIL TAGHIA