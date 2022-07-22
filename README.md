# Trajectory-PMB-EOT-BP

This repository contains the Matlab implementation of the trajectory Poisson multi-Bernoulli filter for multiple extended object tracking using particle belief propagation presented in

https://arxiv.org/abs/2207.10164.

This paper is an extension of the following work

Yuxuan Xia, Karl Granström, Lennart Svensson, Ángel F. García-Femández, and Jason L. Williams (2019, July). Extended Target Poisson multi-Bernoulli Mixture Trackers Based on Sets of Trajectories. In 2019 22nd International Conference on Information Fusion (FUSION) IEEE.

Full text is available at https://arxiv.org/abs/1911.09025.

The MATLAB implementation is adapted from https://github.com/meyer-ucsd/EOT-TSP-21.

The estimation performance of the current set of objects is evaluated using the generalised optimal subpattern-assignment (GOSPA) integrated with the Gaussian Wasserstein distance in

Rahmathullah, Abu Sajana, Ángel F. García-Fernández, and Lennart Svensson. "Generalized optimal sub-pattern assignment metric." 2017 20th International Conference on Information Fusion (Fusion). IEEE, 2017.

Video on GOSPA: https://www.youtube.com/watch?v=M79GTTytvCM

The estimation performance of the set of trajectories is evaluated using the trajectory metric integrated with the Gaussian Wasserstein distance in 

García-Fernández, Ángel F., Abu Sajana Rahmathullah, and Lennart Svensson. "A metric on the space of finite sets of trajectories for evaluation of multi-target tracking algorithms." IEEE Transactions on Signal Processing 68 (2020): 3917-3928.
