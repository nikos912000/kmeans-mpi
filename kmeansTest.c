#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "kmeans.h"
#include "cluster.h"
#include <sys/time.h>
#include <math.h>
#include <mpi.h>


#define max_iterations 50
#define threshold 0.001

void error_message()
{
	char *help = "Error using kmeans: Three arguments required\n"
	"First: number of elements\n"
	"Second: number of attributes (dimensions)\n"
	"Third: numder of clusters\n";

	printf("%s", help);
}

void random_initialization(data_struct *data_in)
{
	int i, j = 0;
	int n = data_in->leading_dim;
	int m = data_in->secondary_dim;
	double *tmp_dataset = data_in->dataset;
	unsigned int *tmp_Index = data_in->members;

	//srand(time(NULL));
	srand(0); 

	for (i = 0; i < m; i++)
	{
		tmp_Index[i] = 0;
		for (j = 0; j < n; j++)
		{
			tmp_dataset[i * n + j] = (double) rand() / RAND_MAX; 
    	}
  	}
}


void initialize_clusters(data_struct *data_in,data_struct *cluster_in)
{
	int i, j, pick = 0;
	int n = cluster_in->leading_dim;
	int m = cluster_in->secondary_dim;
	int Objects = data_in->secondary_dim;
	double *tmp_Centroids = cluster_in->dataset;
	double *tmp_dataset = data_in->dataset;
	unsigned int *tmp_Sizes = data_in->members;

	int step = Objects / m;

	/*randomly pick initial cluster centers*/
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
      		tmp_Centroids[i * n + j] = tmp_dataset[pick * n + j];
    	}
		pick += step; 
	}	
}

void print(data_struct* data2print)
{
	int i, j = 0;
	int n = data2print->leading_dim;
	int m = data2print->secondary_dim;
	double *tmp_dataset = data2print->dataset;

  
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
      		printf("%f ", tmp_dataset[i * n + j]);
    	}
    	printf("\n");
  	}
}


void save(data_struct* data2save, char *filename1, char *filename2)
{
	int i, j = 0;
	FILE *outfile;
	int n = data2save->leading_dim;
	int m = data2save->secondary_dim;
	double *tmp_dataset = data2save->dataset;
	unsigned int *tmp_members = data2save->members;

	printf("Saving data to files: "); 
	printf("%s and %s \n", filename1, filename2); 

	/*===========Save to file 1===========*/
	if((outfile=fopen(filename1, "wb")) == NULL)
	{
    	printf("Can't open output file\n");
  	}

	fwrite(tmp_dataset, sizeof(double), m * n, outfile);

	fclose(outfile);

	/*===========Save to file 2========*/

	if ((outfile = fopen(filename2, "wb")) == NULL)
	{
    	printf("Can't open output file\n");
  	}

	fwrite(tmp_members, sizeof(unsigned int), m, outfile);

	fclose(outfile);
}

void clean(data_struct* data1){

  free(data1->dataset);
  free(data1->members);
}


int main(int argc, char **argv)
{
	int i, processors, rank; 
	struct timeval first, second, lapsed;
	struct timezone tzp;

	

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &processors);
	
	if (argc<4 && rank == 0)
	{
		error_message();
		return 0;
	}

	int numObjects = atoi(argv[1]);
	int numAttributes = atoi(argv[2]);
	int numClusters = atoi(argv[3]);

	char *file1_0 = "centroids.bin";
	char *file1_1 = "ClusterSize.bin";
	char *file2_0 = "dataset.bin";
	char *file2_1 = "Index.bin"; 

	data_struct data_in;
	data_struct clusters;

	/*=======Memory Allocation=========*/
	data_in.leading_dim = numAttributes;
	data_in.secondary_dim = numObjects;
	data_in.dataset = (double*)malloc(numObjects * numAttributes * sizeof(double));
	data_in.members = (unsigned int*)malloc(numObjects * sizeof(unsigned int));

	clusters.leading_dim = numAttributes;
	clusters.secondary_dim = numClusters;
	clusters.dataset = (double*)malloc(numClusters * numAttributes * sizeof(double));
	clusters.members = (unsigned int*)malloc(numClusters * sizeof(unsigned int)); 


	/*=============initialize ==========*/
	random_initialization(&data_in);

	if (rank == 0)
	{
		initialize_clusters(&data_in, &clusters);
 	 }
	
	MPI_Barrier(MPI_COMM_WORLD);

	// Broadcasts a message from the process with rank "root" to all other processes of the communicator 	
	MPI_Bcast(clusters.dataset, numClusters*numAttributes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Create a struct for each process. Each process has its own dataset (a subdata of the initial one)
	data_struct p_data;

	p_data.leading_dim = numAttributes; 

	// The elements associated with each process are n/p
	double n_split =  numObjects / processors;
	double p_objects = ceil(n_split);
	int n_temp = p_objects *  processors;
	if (rank != 0)
	{
		p_data.secondary_dim = p_objects;
		p_data.dataset = (double*)malloc(p_objects * numAttributes * sizeof(double));
		p_data.members = (unsigned int*)malloc(p_objects * sizeof(unsigned int)); 
	}
	else
	{
		p_data.secondary_dim = p_objects + (numObjects - n_temp);
		p_data.dataset = (double*)malloc(p_data.secondary_dim * numAttributes * sizeof(double));
		p_data.members = (unsigned int*)malloc(p_data.secondary_dim * sizeof(unsigned int)); 	
	}
	printf("\nI am process %d and my size is = %d", rank, p_data.secondary_dim);

	// Here, every process creates a sub-dataset of its elements
	for (i = 0; i < p_data.secondary_dim * p_data.leading_dim;i++)
	{ 
		p_data.dataset[i] = data_in.dataset[rank * p_data.secondary_dim * p_data.leading_dim+i]; 
	}

	// Get timestamp (only once!)
	if (rank == 0)
	{
		gettimeofday(&first, &tzp);
	}

	/************************* Clustering *************************/

	/***** Initializations *****/
	int iter, j, k;

	// SumOfDist is the sum of the old iteration 
	// new_SumOfDist is the sum of the current iteration 
	// The difference of them is compared with the threshold
	double SumOfDist = 0, new_SumOfDist=0, temp_SumOfDist=0;
	double* newCentroids;
	int* temp_clusterSize;
	unsigned int*temp_dataMembers;

	// Temporary array for cluster sizes 
	// Used in reduction, as we can't use the same argument for send and receive
	temp_clusterSize = (int*)malloc(numClusters * sizeof(int));

	// Temporary array for data members -in which cluster they belong to- (used in reduction)
	temp_dataMembers = (unsigned int*)malloc(numObjects * sizeof(unsigned int));

	// Our new centroids (used in kmeans_process)
	newCentroids = (double*)malloc(numAttributes * numClusters * sizeof(double));

	// Initialize all cluster sizes to zero
	for (i = 0; i < numClusters; i++)
	{
		temp_clusterSize[i] = 0;
	}

	// Make sure that all processes start clustering at the same time
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0)
	{
		printf("\n\nClastering is going to start!\n");
	}
	/***** Start clustering *****/
	for (iter = 0; iter < max_iterations; iter++)
	{
		new_SumOfDist = 0;
		temp_SumOfDist = 0;

		for (i = 0; i < numClusters * numAttributes; i++)
		{
			newCentroids[i] = 0;
		}

		kmeans_process(&p_data, &clusters, newCentroids, &new_SumOfDist);

		// "Pass" centroids to the struct
		MPI_Allreduce(newCentroids, clusters.dataset, numClusters*numAttributes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		// Each process has its own cluster sizes, so we have to sum them up!
		// Actually, we could do that only once, at the end of clustering
		MPI_Allreduce(clusters.members, temp_clusterSize, numClusters, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);  

		// "Pass" cluster sizes to the struct
		for (i = 0; i < numClusters; i++)
		{
       		clusters.members[i] = temp_clusterSize[i];
		}

		// Update cluster centers
		for (i = 0; i < numClusters; i++)
		{
			for (j = 0; j < numAttributes; j++)
			{
				clusters.dataset[i * numAttributes + j] /= (double) clusters.members[i];
			}
		}
		
		// Calculate total sum from local sums of processes
		MPI_Allreduce(&new_SumOfDist, &temp_SumOfDist, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		// If we have reached threshold, stop clustering
		if (fabs(SumOfDist - temp_SumOfDist) < threshold)
		{
			break;
		}

		// Get the new sum
		SumOfDist = temp_SumOfDist;

		// And print it (once only)
		if (rank == 0)
		{
			printf("\nSum of Distances of iteration %d: %f", iter + 1, SumOfDist);
		}
	}

	// Free some space
	free(newCentroids);
	free(temp_clusterSize);
	
	MPI_Barrier(MPI_COMM_WORLD);

	/* We need some code for saving data members now*/
	
	// Each process saves its own data in a temporary array
	for (i = 0; i < p_data.secondary_dim; i++) 
	{ 
		temp_dataMembers[rank * p_data.secondary_dim + i] = p_data.members[i]; 
	} 

	// Merge the above data
	MPI_Allreduce(temp_dataMembers, data_in.members, numObjects, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

	free(temp_dataMembers);

	// Wait for all processes to reach this point
	MPI_Barrier(MPI_COMM_WORLD);

	// We have finished!
	// Print some stuff and save results in files
	if (rank == 0)
	{
		printf("\n\nFinished after %d iterations\n", iter);
  
		gettimeofday(&second, &tzp);

		if (first.tv_usec > second.tv_usec)
		{
			second.tv_usec += 1000000;
			second.tv_sec--;
		}
  
		lapsed.tv_usec = second.tv_usec - first.tv_usec;
		lapsed.tv_sec = second.tv_sec - first.tv_sec;

		printf("Time elapsed: %d.%06dsec\n\n", (int)lapsed.tv_sec, (int)lapsed.tv_usec); 

		printf("Cluster sizes\n");
		for (i = 0; i < numClusters; i++)
		{
       		printf("Cluster%d: %d\n", i,clusters.members[i]);
		}
		printf("\n");

		/*========save data============*/
		save(&clusters, file1_0, file1_1);
		save(&data_in, file2_0, file2_1);
		printf("\nSaved!\n");
	}

	/*============clean memory===========*/
	clean(&p_data);	
	clean(&data_in);
	clean(&clusters);

	if (rank == 0)
	{
		printf("Program has finished!\n");
	}
	
	MPI_Finalize();
}
