#include "partition-validation.h"

using namespace std;
using std::cout;
using std::cin;
using std::endl;

unsigned stou(char *s){
  return strtoul(s,(char **)NULL,10);
}

  // Call: trade <seed> <Ntries>
int main(int argc,char *argv[]){

  cout << "Version: December 12, 2019." << endl;
  cout << "Command: ";
  cout << argv[0];
  for(int i=1;i<argc; i++)
    cout << " " << argv[i];
  cout << endl;

  // Parse command input
  const string CALL_SYNTAX = "Call: ./partition-validation [-h] [-s <seed>] [-t <distance threshold>] [--validate <validation size>] [--k-fold-crossvalidate k] [--validation-sampling <validation size> <validation samples>] input_partitions.txt output_clustering_txt\n";
  if( argc == 1 ){
    cout << CALL_SYNTAX;
    exit(-1);
  }
  unsigned int seed = 1234;

  string inFileName = "noname";
  string outFileName = "noname";
  int Nskiplines = 0;

  int argNr = 1;
  double distThreshold = 0.2;
  int NvalidationPartitions = 0;
  int NtrainingPartitions = 0;
  int crossvalidateK = 0;
  int validationSamples = 0;
  while(argNr < argc){
    if(to_string(argv[argNr]) == "-h"){
      cout << CALL_SYNTAX;
      cout << "-h: This help" << endl;
      cout << "seed: Any positive integer." << endl; 
      cout << "distance threshold: The max distance from center partition to other partitions in the cluster. Default is 0.2." << endl;
      cout << "--validate <validation size>: The number of partitions at the end of the input that will be used for validation. The first partitions will be used to find clusters. Default is 0 validation partitions." << endl; 
      cout << "--k-fold-crossvalidate k: Perform k-fold cross-validation of all partitions. The training partitions will be used to find clusters and the other ones for validation. Default is 0 folds for no cross-validation." << endl;
      cout << "--validation-sampling <training size> <validation size> <validation samples>: A random set of validation size partitions will be held out from training size partitions in each of validation samples resamplings. Reports the average fraction of validation partitions that belong to clusters. Default is no validation sampling." << endl; 
      cout << "input_partitions.txt: Each column corresponds to a partition and each row corresponds to a node id." << endl;  
      cout << "output_clustering.txt: clusterID partitionID" << endl;  
      exit(-1);
    }
    else if(to_string(argv[argNr]) == "-s"){
      argNr++;
      seed = atoi(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "--skiplines"){
      argNr++;
      Nskiplines = atoi(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-t"){
      argNr++;
      distThreshold = atof(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "--validate"){
      argNr++;
      NvalidationPartitions = atoi(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "--k-fold-crossvalidate"){
      argNr++;
      crossvalidateK = atoi(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "--validation-sampling"){
      argNr++;
      NtrainingPartitions = atoi(argv[argNr]);
      argNr++;
      NvalidationPartitions = atoi(argv[argNr]);
      argNr++;
      validationSamples = atoi(argv[argNr]);
      argNr++;
    }       
    else{

      if(argv[argNr][0] == '-'){
        cout << "Unknown command: " << to_string(argv[argNr]) << endl;
        cout << CALL_SYNTAX;
        exit(-1);
      }

      inFileName = string(argv[argNr]);
      argNr++;
      outFileName = string(argv[argNr]);
      argNr++;
    }

  }

  if(inFileName == "noname"){
    cout << "Missing infile" << endl;
    cout << CALL_SYNTAX;
    exit(-1);
  }
  if(outFileName == "noname"){
    cout << "Missing outfile" << endl;
    cout << CALL_SYNTAX;
    exit(-1);
  }
  
  cout << "Setup:" << endl;
  cout << "-->Using seed: " << seed << endl;
  cout << "-->Will cluster partitions such that no partition is farther away from its center than: " << distThreshold << endl; 
  cout << "-->Will read partitions from file: " << inFileName << endl;
  if(crossvalidateK > 0){
    cout << "-->Will perform " << crossvalidateK << "-fold cross-valiation." << endl;
    NvalidationPartitions = 0;
    validationSamples = 0;
  }
  else if(validationSamples > 0){
    cout << "-->Will perform validation sampling of " << NvalidationPartitions << " partitions given " << NtrainingPartitions << " training partitions " << validationSamples << " times." << endl;
    crossvalidateK = 0;
  }   
  else if(NvalidationPartitions > 0){
    cout << "-->Will use the last " << NvalidationPartitions << " partitions for validation." << endl;
    crossvalidateK = 0;
    validationSamples = 0;
  }
 
  cout << "-->Will write clusters to file: " << outFileName << endl;
  cout << "-->Will use number of threads: " <<  omp_get_max_threads() << endl;

  Partitions partitions(inFileName,outFileName,distThreshold,NtrainingPartitions,NvalidationPartitions,crossvalidateK,validationSamples,seed);

  if(crossvalidateK == 0){
    if(validationSamples == 0){
      partitions.clusterPartitions(0);
      if(NvalidationPartitions > 0){
        partitions.validatePartitions(0,"");
        cout << "-->Fraction of validation partitions that fits in a cluster: " << 1.0*partitions.NtotValidated/partitions.NtotTested << endl;
      }
      partitions.printClusters();
    }
    else{
      for(int sample = 0; sample < validationSamples; sample++){
        partitions.clusterPartitions(0);
        partitions.validatePartitions(0,"");
        cout << "-->Fraction of validation partitions that fits in a cluster after " << sample+1 << " samples: " << 1.0*partitions.NtotValidated/partitions.NtotTested << endl;
      }
    }
  }
  else{
    for(int fold = 0;fold<crossvalidateK;fold++){
      cout << endl << "Fold " << fold+1 << "/" << crossvalidateK << endl;
      partitions.clusterPartitions(fold);
      partitions.validatePartitions(fold,"_" + to_string(fold+1));
      cout << "-->Fraction of validation partitions that fits in a cluster after " << fold+1 << " folds: " << 1.0*partitions.NtotValidated/partitions.NtotTested << endl;
    }

  }

}