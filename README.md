Parallel KMeans using MPI
=========================

A parallel implementation of KMeans, in C, using MPI.

2nd Course Assignment for Parallel and Distributed Computing Systems (2013).

How to use
----------

1. Install mpich2: `sudo apt-get install libcr-dev mpich2 mpich2-doc`
2. Run `make` command in a unix-based system
3. Run the executable using `mpirun -np N`, where N is the number of processes.

#####Arguments#####

*number of points, dimensions, clusters*  

Example: 

For 3 processes, 100 points, 3 dimensions and 3 clusters:

`mpirun -np 3 ./kmeansTest 100 3 3`

Output
------
Info messages in stdout. Results are saved to files (centroids, ClusterSize, dataset and Index).