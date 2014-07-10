#include "kmeans.h"
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

double euclidean_distance(double *v1, double *v2, int length)
{
	int i = 0;
	double dist = 0;

	for (i = 0; i < length; i++)
	{
		dist += (v1[i] - v2[i]) * (v1[i] - v2[i]); 
	 }

	return(dist);
}


void kmeans_process(data_struct *data_in, data_struct *clusters, double *newCentroids, double* SumOfDist)
{
	int i, j, k;
	double tmp_dist = 0;
	int tmp_index = 0;
	double min_dist = 0;
	double *dataset = data_in->dataset;
	double *centroids = clusters->dataset;
	unsigned int *Index = data_in->members;
	unsigned int *cluster_size = clusters->members;

	for (i = 0; i < clusters->secondary_dim; i++)
	{
		cluster_size[i] = 0;
	}

	for (i = 0; i < data_in->secondary_dim; i++)
	{
		tmp_dist = 0;
		tmp_index = 0;
		min_dist = FLT_MAX;
		
		/*find nearest center*/
		for (k = 0; k < clusters->secondary_dim; k++)
		{
			tmp_dist = euclidean_distance(dataset + i * data_in->leading_dim, centroids + k * clusters->leading_dim, data_in->leading_dim);
     
			if (tmp_dist<min_dist)
			{
				min_dist = tmp_dist;
				tmp_index = k;
			}
		}
	   
		Index[i] = tmp_index;
		SumOfDist[0] += min_dist;
		cluster_size[tmp_index]++;
		for (j = 0; j < data_in->leading_dim; j++)
		{
			newCentroids[tmp_index * clusters->leading_dim + j] += dataset[i * data_in->leading_dim + j]; 
		}
	}
}