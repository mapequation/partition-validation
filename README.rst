Hierarchical partition validation
=================================

Identifies partition clusters with partitions within a given distance from their centers. Measures distances by a weighted Jaccard index over possibly hierarchical partitions. Optionally validates how many partitions in a random or selected hold-out subset that fall within the clusters.

Preserves the partition order in the input file such that if they are given in decreasing quality, the center of each cluster will have the highest quality within the cluster.

Complexity is linear in the number of partitions. 

Getting started
---------------

In a terminal with the GNU Compiler Collection installed,
just run ``make`` in the current directory to compile the
code with the included Makefile.


Run the binary ``partition-validation`` using::

    ./partition-validation [-h] [-s <seed>] [-t <distance threshold>] [--validate N] [--k-fold-crossvalidate <k>] input_partitions.txt output_clustering_txt

For a list of options, run::

    ./partition-validation -h
    -h: This help.  
    seed: Any positive integer.  
    distance threshold: The max distance from center partition to other partitions in the cluster. Default is 0.2.  
    --validate N: The number of partitions N at the end of the input that will be used for validation. The first partitions will be used to find clusters. Default is 0 validation partitions.  
    --k-fold-crossvalidate k: Perform k-fold cross-validation of all partitions. The training partitions will be used to find clusters and the other ones for validation. Default is 0 folds for no cross-validation.  
    --validation-sampling <training size> <validation size> <validation samples>: A random set of validation size partitions will be held out from training size partitions in each of validation samples resamplings. Reports the average fraction of validation partitions that belong to clusters. Default is no validation sampling.  
    input_partitions.txt: Each column corresponds to a partition and each row corresponds to a node id.  
    output_clustering.txt: clusterID partitionID.  
    validation_output_clustering.txt: bool (1 if validation partition fits in a cluster and 0 if not, printed in the order of the validation partitions)  
  
Examples
~~~~~~~~

To cluster 10 hierarchical partitions in input file:: 

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

with distance threshold 0.2, run::

    ./partition-validation -s 123 -t 0.2 input_partitions.txt output_clustering.txt 

and generate output files::

    output_clustering.txt   
	# Clustered 10 partitions into 2 clusters.
	# ClusterId PartitionId
	1 1
	1 2
	1 3
	1 4
	1 5
	2 6
	2 7
	2 8
	2 9
	2 10

and::

	output_clustering_distances.txt
	# ClusterId1 ClusterId2 Distance
	1 2 0.398148148148148     

To validate the same hierarchical partitions with 5 training and 5 validation partitions sampled 10 times, run::
	
	./partition-validation -s 123 -t 0.2 --validation-sampling 5 5 10 input_partitions.txt partitions_clustering.txt

and generate output file::

	partitions_clustering_validation.txt
	# number_of_validated_out_of_5 number_of_partition_clusters
	2 5
	4 4
	3 3
	2 2
	3 3
	2 5
	3 3
	3 3
	5 2
	5 2


Author
------

Martin Rosvall

For contact information, see: http://www.mapequation.org/about  
