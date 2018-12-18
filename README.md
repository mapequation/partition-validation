# Hierarchical partition validation

Identifies partition clusters with partitions within a given distance from their centers. Measures distances by a weighted Jaccard index over hierarchical possibly partitions. Optionally validates how many partitions in a random or selected hold-out subset that fall within the clusters.

Preserves the partition order in the input file such that if they are given in decreasing quality, the center of each cluster will have the highest quality within the cluster.

Complexity is linear in the number of partitions. 

## Author:

Martin Rosvall

For contact information, see http://www.mapequation.org/about

## Getting started:

In a terminal with the GNU Compiler Collection installed,
just run 'make' in the current directory to compile the
code with the included Makefile.


Call: ./partition-validation [-h] [-s \<seed\>] [-t \<distance threshold\>] [--skiplines N] [--validate N] [--k-fold-crossvalidate \<k\>] input_partitions.txt output_clustering_txt  

-h: This help.  
seed: Any positive integer.  
distance threshold: The max distance from center partition to other partitions in the cluster. Default is 0.2.  
--validate N: The number of partitions N at the end of the input that will be used for validation. The first partitions will be used to find clusters. Default is 0 validation partitions. 
--k-fold-crossvalidate k: Perform k-fold cross-validation of all partitions. The training partitions will be used to find clusters and the other ones for validation. Default is 0 folds for no cross-validation.
--validation-sampling \<validation size\> \<validation samples\>: A random set of validation size partitions will be held out in each of validation samples resamplings. Reports the average fraction of validation partitions that belong to clusters. Default is no validation sampling.    
input_partitions.txt: Each column corresponds to a partition and each row corresponds to a node id.  
output_clustering.txt: clusterID partitionID.
validation_output_clustering.txt: bool (1 if validation partition fits in a cluster and 0 if not, printed in the order of the validation partitions)  
  
Example:

./partition-validation -s 123 -t 0.2 input_partitions.txt output_clustering.txt    

input_multilevel_partitions.txt  
1:1 1:1 1:1 1:1 1:1 1:1 1:1 1:1 1:1 1:1  
1:1 1:1 1:1 1:1 1:2 1:1 1:1 1:1 1:2 1:1  
1:2 1:2 1:2 1:2 1:3 1:2 1:1 1:2 1:2 1:2  
1:2 1:3 1:3 1:2 1:4 1:2 1:2 1:2 1:2 1:2  
2:1 2:1 2:1 2:1 2:1 1:3 1:2 1:2 1:3 1:3  
2:2 2:2 2:2 2:2 2:1 2:1 2:1 2:1 2:1 2:1  
3:1 3:1 3:1 3:1 3:1 2:1 2:1 2:1 2:1 2:1  
3:1 3:1 3:1 3:2 3:2 2:2 2:1 2:2 2:1 2:1  
3:2 3:2 3:2 3:2 3:3 2:2 2:2 2:2 2:2 2:2  

output_clustering.txt   
\# Clustered 10 partitions into 2 clusters.  
\# ClusterId PartitionId  
\# Cluster 1: 5 partitions.  
1 1  
1 2  
1 3  
1 4  
1 5  
\# Cluster 2: 5 partitions.  
2 6  
2 7  
2 8  
2 9  
2 10   
