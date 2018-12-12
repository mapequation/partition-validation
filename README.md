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
distance threshold: The max distance between two partitions in any cluster. Default is 0.2.  
--skiplines N: Skip N lines in input_partitions.txt before reading data.  
--validate N: The number of partitions N at the end of the input that will be used for validation. The first partitions will be used to find clusters. Default is 0 validation partitions. 
--k-fold-crossvalidate k: Perform k-fold cross-validation of all partitions. The training partitions will be used to find clusters and the other ones for validation. Default is 0 folds for no cross-validation.    
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
\# Clustered 10 partitions into 2 clusters with maximum distance 0.466666666666667, average maximum distance  0.344444444444444,  and maximum cluster size 5  
\# ClusterId PartitionId  
\# Cluster 1: 5 nodes with max distance 0.466666666666667  
1 1  
1 2  
1 3  
1 4  
1 5  
\# Cluster 2: 5 nodes with max distance 0.222222222222222  
2 6  
2 7  
2 8  
2 9  
2 10   
