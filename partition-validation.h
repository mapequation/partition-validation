#include <cstdio>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <random>
#include <functional>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#ifdef _OPENMP
#include <omp.h>
#include <stdio.h>
#else
  #define omp_get_thread_num() 0
  #define omp_get_max_threads() 1
#endif
using namespace std;
#include <limits>
const double epsilon = 1.0e-15;
const double bignum = 1.0;
const double threshold = 1.0e-10;
const unsigned int max_unsigned_int_size = std::numeric_limits<unsigned int>::max();

// ofstream with higher precision to avoid truncation errors
struct my_ofstream : ofstream {
  explicit my_ofstream(streamsize prec = 15)
  {
    this->precision(prec);
  }
};


template <class T>
inline string to_string (const T& t){
	stringstream ss;
	ss << t;
	return ss.str();
}

vector<string> tokenize(const string& str,string& delimiters)
{

  vector<string> tokens;

  // skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);

  // find first "non-delimiter".
  string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos){

    // found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));

    // skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);

    // find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);

  }

  return tokens;
}

struct pairhash {
public:
  template <typename T, typename U>
  size_t operator()(const pair<T, U> &x) const
  {
    return x.first*31 + x.second;
  }
};

unsigned int convert_div_to_uintdiv(double div)
{ 
  return static_cast<unsigned int>(max_unsigned_int_size*div);
}

double convert_uintdiv_to_div(unsigned int uintdiv)
{ 
  return static_cast<double>(1.0*uintdiv/max_unsigned_int_size);
}

class Partition{
public:
	Partition();
	Partition(int partitionid, int Nnodes);
	int partitionId;
	vector<vector<int> > assignments;
	unordered_map<int,int> clusterSizes;
	Partition *minCenterPartition;
};

Partition::Partition(){
};

Partition::Partition(int partitionid, int Nnodes){
	partitionId = partitionid;
	assignments = vector<vector<int> >(Nnodes);
}

typedef vector< vector<Partition *> > Clusters;



class Partitions{
private:

	// double wpJaccardDist(int partitionId1, int partitionId2);
	double wpJaccardDist(Partition *partition1, Partition *partition2);
	
	vector<Partition> partitions;
	int crossvalidateK = 0;

	int Nthreads = 1;

	int randInt(int from, int to);
	double randDouble(double to);

	string inFileName;
	string outFileName;
	double distThreshold = 0.2;
	vector<mt19937> mtRands;
	ifstream ifs;
  	string line;
	Clusters clusters;

public:
	Partitions(string inFileName,string outFileName,double distThreshold,int NvalidationPartitions,int crossvalidateK,int seed); 
	void readPartitionsFile();
	void clusterPartitions(int fold);
	void printClusters();
	void validatePartitions(int fold,string filesuffix);

	int Nnodes = 0;
	int Npartitions = 0;
	int Nclusters = 0;
	int NtrainingPartitions = 0;
	int NvalidationPartitions = 0;

	int NtotValidated = 0;
	int NtotTested = 0;
};

Partitions::Partitions(string inFileName,string outFileName,double distThreshold,int NvalidationPartitions,int crossvalidateK,int seed){
	this->distThreshold = distThreshold;
	this->NvalidationPartitions = NvalidationPartitions;
	this->crossvalidateK = crossvalidateK;
	this->inFileName = inFileName;
	this->outFileName = outFileName;

	Nthreads = max(1, omp_get_max_threads());
	for(int i = 0; i < Nthreads; i++)
  	  mtRands.push_back(mt19937(seed+1));
  
	// Open state network for building state node to physical node mapping
	ifs.open(inFileName.c_str());
	if(!ifs){
		cout << "failed to open \"" << inFileName << "\" exiting..." << endl;
		exit(-1);
	}
	readPartitionsFile();
	ifs.close();

}

int Partitions::randInt(int from, int to){

	uniform_int_distribution<int> rInt(from,to);
	return rInt(mtRands[omp_get_thread_num()]);

}

double Partitions::randDouble(double to){

	uniform_real_distribution<double> rDouble(0.0,to);
	return rDouble(mtRands[omp_get_thread_num()]);

}

double Partitions::wpJaccardDist(Partition *partition1, Partition *partition2){

	int partitionId1Size = partition1->clusterSizes.size();
	int partitionId2Size = partition2->clusterSizes.size();
	vector<double> maxClusterSimilarityPartitionId1(partitionId1Size,0.0);
	vector<double> maxClusterSimilarityPartitionId2(partitionId2Size,0.0);

	unordered_map<pair<int,int>,int,pairhash> jointM;
	for(int k=0;k<Nnodes;k++){
		for(vector<int>::iterator it_id1 = partition1->assignments[k].begin(); it_id1 != partition1->assignments[k].end(); it_id1++){
	  	for(vector<int>::iterator it_id2 = partition2->assignments[k].begin(); it_id2 != partition2->assignments[k].end(); it_id2++)
	  		jointM[make_pair((*it_id1),(*it_id2))]++;
	  }
	}
	
	for(unordered_map<pair<int,int>,int,pairhash>::iterator it = jointM.begin(); it != jointM.end(); it++){

	  int Ncommon = it->second;
	  int Ntotal = partition1->clusterSizes[it->first.first] + partition2->clusterSizes[it->first.second] - Ncommon;
	  double clusterSim = 1.0*Ncommon/Ntotal;
	  maxClusterSimilarityPartitionId1[it->first.first] = max(clusterSim,maxClusterSimilarityPartitionId1[it->first.first]);
	  maxClusterSimilarityPartitionId2[it->first.second] = max(clusterSim,maxClusterSimilarityPartitionId2[it->first.second]);

	}

	int NassignmentsId1 = 0;
	double simId1 = 0.0;
	for(int i=0;i<partitionId1Size;i++){
		simId1 += maxClusterSimilarityPartitionId1[i]*partition1->clusterSizes[i];
		NassignmentsId1 += partition1->clusterSizes[i];
	}
	simId1 /= 1.0*NassignmentsId1;

	int NassignmentsId2 = 0;
	double simId2 = 0.0;
	for(int i=0;i<partitionId2Size;i++){
		simId2 += maxClusterSimilarityPartitionId2[i]*partition2->clusterSizes[i];
		NassignmentsId2 += partition2->clusterSizes[i];
	}
	simId2 /= 1.0*NassignmentsId2;

	return 1.0 - 0.5*simId1 - 0.5*simId2;

}

void Partitions::clusterPartitions(int fold){

	#ifdef _OPENMP
	// Initiate locks to keep track of best solutions
	omp_lock_t lock;
	omp_init_lock(&lock);
	#endif

	cout << "Clustering partitions:" << endl;


	vector<Partition*> partitionPtrs = vector<Partition*>(NtrainingPartitions);
	for(int i=0;i<NtrainingPartitions;i++){
		partitionPtrs[i] = &partitions[i + (i >= (Npartitions-(fold+1)*NvalidationPartitions) ? NvalidationPartitions : 0) ];
	}
	clusters.push_back({partitionPtrs[0]});
	Nclusters++;
	// vector<Partition *> clusteredPartitionPtrs{partitionPtrs[0]};
	// fastClusters.push_back(clusteredPartitionPtrs);
	for(int i=1;i<NtrainingPartitions;i++){
		bool fitsInCluster = false;
		for(int j=0;j<Nclusters;j++){
			if(wpJaccardDist(partitionPtrs[i],clusters[j][0]) < distThreshold){
				fitsInCluster = true;
				clusters[j].push_back(partitionPtrs[i]);
				break;
			}
		}
		if(!fitsInCluster){
			clusters.push_back({partitionPtrs[i]});
			Nclusters++;
		}
	}


}

void Partitions::validatePartitions(int fold,string filesuffix){

	int Nvalidated = 0;

	vector<Partition*> validationPartitionPtrs = vector<Partition*>(NvalidationPartitions);
	for(int i=0;i<NvalidationPartitions;i++){
		validationPartitionPtrs[i] = &partitions[Npartitions-(fold+1)*NvalidationPartitions+i];
	}
	cout << "-->Number of validation partitions out of " << NvalidationPartitions << " (" << validationPartitionPtrs[0]->partitionId+1 << "-" << validationPartitionPtrs[NvalidationPartitions-1]->partitionId+1 << ") that fits in a cluster is..." << flush;
	vector<int> validatedPartitions(NvalidationPartitions,0);

	#pragma omp parallel for
	for(int i=0;i<NvalidationPartitions;i++){
		
		for(int j=0;j<Nclusters;j++){
			if(wpJaccardDist(validationPartitionPtrs[i],clusters[j][0]) < distThreshold){
				#pragma omp atomic
				Nvalidated++;
				validatedPartitions[i] += 1;				
				break;
			}
		}
	}
	cout << Nvalidated << ". " << flush;	
	NtotValidated += Nvalidated;
	NtotTested += NvalidationPartitions;

	string validationOutFileName = outFileName;
	size_t period_pos = validationOutFileName.find_last_of(".");
	if(period_pos ==  string::npos)
		validationOutFileName += "_validation" + filesuffix;
	else
		validationOutFileName.insert(period_pos,"_validation" + filesuffix);
	cout << "Writing validation results to " << validationOutFileName << endl;
	my_ofstream ofs;
	ofs.open(validationOutFileName.c_str());
	for(int i=0;i<NvalidationPartitions;i++){
		ofs << validationPartitionPtrs[i]->partitionId+1 << " " << validatedPartitions[i] << endl;
	}
	ofs.close();

}

void Partitions::readPartitionsFile(){

  cout << "Reading partitions file " << flush;  

  string line;  
  string buf;
  int lineNr = 0;

  // Count number of nodes and boot partitions
  getline(ifs,line);
  lineNr++;
  Nnodes++; // First line corresponds to first node

  istringstream read(line);
  while(read >> buf)
      Npartitions++;

  if(crossvalidateK > 0)
  	NvalidationPartitions = Npartitions/crossvalidateK;

  if(Npartitions - NvalidationPartitions > 0){
  	NtrainingPartitions = Npartitions - NvalidationPartitions;
  }
  else{
  	cout << "--Not enough partitions for validation. Will not validate.--" << endl;
  }

  cout << "with " << Npartitions << " partitions with " << NtrainingPartitions << " partitions for clustering and " << NvalidationPartitions << " partitions for validation " << flush;
  if(crossvalidateK > 0)
  	cout << "in each of the " << crossvalidateK << " folds " << flush;

  // Count remaining nodes
  while(getline(ifs,line))
    Nnodes++;
  if(ifs.bad())
    cout << "Error while reading file" << endl;

  cout << "of " << Nnodes << " nodes..." << flush;

  partitions = vector<Partition>(Npartitions);
  for(int i=0;i<Npartitions;i++)
  	partitions[i] = Partition(i,Nnodes);

  // Restart from beginning of file
  vector<map<string,int> > partitionsAssignmentId(Npartitions);
  vector<int> partitionsAssignmentIds(Npartitions,0);
  string delim = ":";
  ifs.clear();
  ifs.seekg(0, ios::beg);

  // Read partitions data    
  int nodeNr = 0;
  lineNr = 0;
  while(getline(ifs,line)){
    lineNr++;
    istringstream read(line);
    int i = 0;
    while(read >> buf){
    	string assignment = buf.c_str();
    	vector<string> assignments = tokenize(assignment,delim);
    	string assignmentKey = "";
    	for(vector<string>::iterator it = assignments.begin(); it != assignments.end(); it++){
    		assignmentKey += (*it) + ":";
    		map<string,int>::iterator assignmentId_it = partitionsAssignmentId[i].find(assignmentKey);
    		int assignmentId = partitionsAssignmentIds[i];
    		if(assignmentId_it != partitionsAssignmentId[i].end()){			assignmentId = assignmentId_it->second;
    		}
    		else{
    			partitionsAssignmentId[i][assignmentKey] = partitionsAssignmentIds[i];
    			partitionsAssignmentIds[i]++;
    		}
    		partitions[i].assignments[nodeNr].push_back(assignmentId); 
      	partitions[i].clusterSizes[assignmentId]++;
    	}
      i++;
    }
    nodeNr++;
  }

  cout << "done!" << endl;

}




void Partitions::printClusters(){

	cout << "-->Writing clustering results..." << flush;

  	my_ofstream ofs;
	ofs.open(outFileName.c_str());

	int i = 1;
	ofs << "# Clustered " << NtrainingPartitions << " partitions into " << clusters.size() << " clusters." << endl;
	ofs << "# ClusterId PartitionId" << endl;
	for(Clusters::iterator cluster_it = clusters.begin(); cluster_it != clusters.end(); cluster_it++){
		ofs << "# Cluster " << i << ": " << cluster_it->size() << " partitions." << endl;
		for(vector<Partition *>::iterator partition_it = cluster_it->begin(); partition_it != cluster_it->end(); partition_it++)
			ofs << i << " " << (*partition_it)->partitionId+1 << endl;
		
		i++;
	}

	ofs.close();

	cout << "done!" << endl;

}

