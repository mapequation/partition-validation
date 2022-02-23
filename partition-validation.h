#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include <random>
#include <functional>
#include <algorithm>
#include <numeric>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#include <stdio.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif


using namespace std;

// ofstream with higher precision to avoid truncation errors
struct my_ofstream : ofstream {
    explicit my_ofstream(streamsize prec = 15) {
        this->precision(prec);
    }
};

template<typename T>
inline string to_string(const T &t) {
    stringstream ss;
    ss << t;
    return ss.str();
}

vector<string> tokenize(const string &str, char delimiter) {
    vector<string> tokens;

    // skip delimiter at beginning.
    auto lastPos = str.find_first_not_of(delimiter, 0);

    // find first "non-delimiter".
    auto pos = str.find_first_of(delimiter, lastPos);

    while (string::npos != pos || string::npos != lastPos) {
        // found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));

        // skip delimiter.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiter, pos);

        // find next "non-delimiter"
        pos = str.find_first_of(delimiter, lastPos);
    }

    return tokens;
}

struct pairhash {
public:
    template<typename T, typename U>
    size_t operator()(const pair<T, U> &x) const {
        return x.first * 31 + x.second;
    }
};

class Partition {
public:
    Partition() = default;

    Partition(int id, int nNodes) : partitionId(id) {
        assignments = vector<vector<int> >(nNodes);
    }

    int partitionId{};
    vector<vector<int> > assignments;
    unordered_map<int, int> clusterSizes;
};


class Partitions {
    double wpJaccardDist(Partition *partition1, Partition *partition2);

    unordered_map<pair<int, int>, double, pairhash> wpJaccardDistLookUp;

    using Clusters = vector<vector<Partition *>>;
    Clusters clusters;

    vector<Partition> partitions;

    int crossValidateK = 0;
    int nThreads = 1;
    string outFileName;
    double distThreshold = 0.2;
    vector<mt19937> mtRands;
    ifstream ifs;
    string validationSampleOutFileName;
    my_ofstream validationOfs;

public:
    int nNodes = 0;
    int nPartitions = 0;
    int nClusters = 0;
    int nTrainingPartitions = 0;
    int nValidationPartitions = 0;
    int validationSamples = 0;
    int nTotValidated = 0;
    int nTotTested = 0;

    Partitions(const string &inFileName,
               const string &outFileName,
               double distThreshold,
               int nTrainingPartitions,
               int nValidationPartitions,
               int crossValidateK,
               int validationSamples,
               int seed)
            : crossValidateK(crossValidateK),
              outFileName(outFileName),
              distThreshold(distThreshold),
              nTrainingPartitions(nTrainingPartitions),
              nValidationPartitions(nValidationPartitions),
              validationSamples(validationSamples) {

        nThreads = max(1, omp_get_max_threads());

        for (int i = 0; i < nThreads; i++) {
            mtRands.emplace_back(seed + 1);
        }

        ifs.open(inFileName.c_str());

        if (!ifs) {
            cout << "Failed to open \"" << inFileName << "\" exiting..." << endl;
            exit(-1);
        }

        readPartitionsFile();

        if (crossValidateK > 0) {
            // Randomize ordered partitions
            shuffle(begin(partitions), end(partitions), mtRands[0]);
        }

        if (validationSamples > 0) {
            validationSampleOutFileName = outFileName;
            auto periodPos = validationSampleOutFileName.find_last_of('.');

            if (periodPos == string::npos) {
                validationSampleOutFileName += "_validation";
            } else {
                validationSampleOutFileName.insert(periodPos, "_validation");
            }

            validationOfs.open(validationSampleOutFileName.c_str());
            validationOfs << "# number_of_validated_out_of_" << nValidationPartitions
                          << " number_of_partition_clusters" << endl;
        }

        ifs.close();
    }

    void readPartitionsFile();

    void clusterPartitions(int fold);

    void printClusters();

    void validatePartitions(int fold);
};

double Partitions::wpJaccardDist(Partition *partition1, Partition *partition2) {
    auto partition1Id = partition1->partitionId;
    auto partition2Id = partition2->partitionId;

    if (partition1Id > partition2Id)
        swap(partition1Id, partition2Id);

    auto lookup = wpJaccardDistLookUp.find(make_pair(partition1Id, partition2Id));

    if (lookup != end(wpJaccardDistLookUp))
        return lookup->second;

    auto partition1Size = partition1->clusterSizes.size();
    auto partition2Size = partition2->clusterSizes.size();

    vector<double> maxClusterSimilarityPartition1(partition1Size, 0.0);
    vector<double> maxClusterSimilarityPartition2(partition2Size, 0.0);

    unordered_map<pair<int, int>, int, pairhash> jointM;

    for (int k = 0; k < nNodes; k++) {
        for (auto id1: partition1->assignments[k]) {
            for (auto &id2: partition2->assignments[k]) {
                jointM[make_pair(id1, id2)]++;
            }
        }
    }

    for (auto &[ids, nCommon]: jointM) {
        int nTotal = partition1->clusterSizes[ids.first] + partition2->clusterSizes[ids.second] - nCommon;
        double clusterSim = 1.0 * nCommon / nTotal;
        maxClusterSimilarityPartition1[ids.first] = max(clusterSim, maxClusterSimilarityPartition1[ids.first]);
        maxClusterSimilarityPartition2[ids.second] = max(clusterSim, maxClusterSimilarityPartition2[ids.second]);
    }

    int nAssignments1 = 0;
    double sim1 = 0.0;

    for (auto i = 0; i < static_cast<int>(partition1Size); i++) {
        sim1 += maxClusterSimilarityPartition1[i] * partition1->clusterSizes[i];
        nAssignments1 += partition1->clusterSizes[i];
    }

    sim1 /= 1.0 * nAssignments1;

    int nAssignments2 = 0;
    double sim2 = 0.0;

    for (auto i = 0; i < static_cast<int>(partition2Size); i++) {
        sim2 += maxClusterSimilarityPartition2[i] * partition2->clusterSizes[i];
        nAssignments2 += partition2->clusterSizes[i];
    }

    sim2 /= 1.0 * nAssignments2;

    double dist = 1.0 - 0.5 * sim1 - 0.5 * sim2;

    wpJaccardDistLookUp[make_pair(partition1Id, partition2Id)] = dist;

    return dist;
}

void Partitions::clusterPartitions(int fold) {
#ifdef _OPENMP
    // Initiate locks to keep track of best solutions
    omp_lock_t lock;
    omp_init_lock(&lock);
#endif

    cout << "Clustering " << nTrainingPartitions << " partitions:" << endl;

    auto partitionPtrs = vector<Partition *>(nTrainingPartitions);

    if (validationSamples > 0) {
        // Shuffle partitions for new training and validation set
        shuffle(begin(partitions), end(partitions), mtRands[0]);
    }

    for (auto i = 0; i < nTrainingPartitions; i++) {
        auto limit = nPartitions - (fold + 1) * nValidationPartitions;
        auto offset = i >= limit ? nValidationPartitions : 0;
        partitionPtrs[i] = &partitions[i + offset];
    }

    if (crossValidateK > 0 || validationSamples > 0) {
        // Order randomized partitions
        sort(begin(partitionPtrs), end(partitionPtrs),
             [](Partition *a, Partition *b) { return a->partitionId < b->partitionId; });
    }

    clusters.clear();
    clusters.push_back({partitionPtrs[0]});
    nClusters++;

    for (int i = 1; i < nTrainingPartitions; i++) {
        auto fitsInCluster = false;
        for (int j = 0; j < nClusters; j++) {
            if (wpJaccardDist(partitionPtrs[i], clusters[j][0]) < distThreshold) {
                fitsInCluster = true;
                clusters[j].push_back(partitionPtrs[i]);
                break;
            }
        }
        if (!fitsInCluster) {
            clusters.push_back({partitionPtrs[i]});
            nClusters++;
        }
    }
}

void Partitions::validatePartitions(int fold) {
    int nValidated = 0;

    auto validationPartitionPtrs = vector<Partition *>(nValidationPartitions);
    for (int i = 0; i < nValidationPartitions; i++) {
        validationPartitionPtrs[i] = &partitions[nPartitions - (fold + 1) * nValidationPartitions + i];
    }

    cout << "--> Number of validation partitions out of " << nValidationPartitions
         << " that fits in one of " << nClusters << " clusters is..." << flush;

    vector<int> validatedPartitions(nValidationPartitions, 0);

#pragma omp parallel for
    for (auto i = 0; i < nValidationPartitions; i++) {
        for (auto j = 0; j < nClusters; j++) {
            if (wpJaccardDist(validationPartitionPtrs[i], clusters[j][0]) < distThreshold) {
#pragma omp atomic
                nValidated++;
                validatedPartitions[i] += 1;
                break;
            }
        }
    }

    cout << nValidated << ". " << endl;
    nTotValidated += nValidated;
    nTotTested += nValidationPartitions;

    if (validationSamples == 0) {
        auto validationOutFileName = outFileName;
        auto periodPos = validationOutFileName.find_last_of('.');
        if (periodPos == string::npos) {
            validationOutFileName += "_validation";
        } else {
            validationOutFileName.insert(periodPos, "_validation");
        }

        cout << "Writing validation results to " << validationOutFileName << endl;

        my_ofstream ofs;
        ofs.open(validationOutFileName.c_str());

        for (int i = 0; i < nValidationPartitions; i++) {
            ofs << validationPartitionPtrs[i]->partitionId + 1 << " " << validatedPartitions[i] << endl;
        }

        ofs.close();
    } else {
        validationOfs << nValidated << " " << nClusters << endl;
    }
}

void Partitions::readPartitionsFile() {
    cout << "Reading partitions file " << flush;

    string line;
    string buf;

    // Count number of nodes and boot partitions
    getline(ifs, line);
    nNodes++; // First line corresponds to first node

    {
        istringstream read(line);
        while (read >> buf) nPartitions++;
    }

    if (crossValidateK > 0) {
        nValidationPartitions = nPartitions / crossValidateK;
    }

    if (nTrainingPartitions == 0) {
        if (nPartitions - nValidationPartitions > 0) {
            nTrainingPartitions = nPartitions - nValidationPartitions;
        } else {
            cout << "-- Not enough partitions for validation. Will not validate." << endl;
        }
    } else if (nPartitions - nValidationPartitions - nTrainingPartitions < 0) {
        cout << "-- Not enough partitions for training and validation." << endl;
    }

    cout << "with " << nPartitions << " partitions..." << flush;

    // Count remaining nodes
    while (getline(ifs, line)) {
        nNodes++;
    }

    if (ifs.bad()) {
        cout << "Error while reading file" << endl;
    }

    cout << "of " << nNodes << " nodes..." << flush;

    partitions = vector<Partition>(nPartitions);

    for (int i = 0; i < nPartitions; i++) {
        partitions[i] = Partition(i, nNodes);
    }

    // Restart from beginning of file
    vector<map<string, int> > partitionsAssignmentId(nPartitions);
    vector<int> partitionsAssignmentIds(nPartitions, 0);
    ifs.clear();
    ifs.seekg(0, ios::beg);

    // Read partitions data
    int nodeNr = 0;

    while (getline(ifs, line)) {
        istringstream read(line);
        int i = 0;

        while (read >> buf) {
            vector<string> assignments = tokenize(buf, ':');
            string assignmentKey;

            for (auto &assignment: assignments) {
                assignmentKey += assignment + ":";

                auto lookup = partitionsAssignmentId[i].find(assignmentKey);
                int assignmentId = partitionsAssignmentIds[i];

                if (lookup != end(partitionsAssignmentId[i])) {
                    assignmentId = lookup->second;
                } else {
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


void Partitions::printClusters() {
    cout << "--> Writing " << nClusters << " clusters..." << flush;

    my_ofstream ofs;
    ofs.open(outFileName.c_str());

    int i = 1;
    ofs << "# Clustered " << nTrainingPartitions << " partitions into " << clusters.size() << " clusters." << endl;
    ofs << "ClusterId PartitionId" << endl;
    for (auto &cluster: clusters) {
        for (auto &partition: cluster) {
            ofs << i << " " << partition->partitionId + 1 << endl;
        }
        i++;
    }

    ofs.close();

    string distanceOutFileName = outFileName;
    size_t periodPos = distanceOutFileName.find_last_of('.');

    if (periodPos == string::npos) {
        distanceOutFileName += "_distances";
    } else {
        distanceOutFileName.insert(periodPos, "_distances");
    }

    ofs.open(distanceOutFileName.c_str());
    i = 1;

    ofs << "# Clustered " << nTrainingPartitions << " partitions into " << clusters.size() << " clusters." << endl;
    ofs << "ClusterId1 ClusterId2 Distance" << endl;

    for (auto cluster1 = begin(clusters); cluster1 != end(clusters); cluster1++) {
        auto j = i + 1;
        for (auto cluster2 = next(cluster1); cluster2 != end(clusters); cluster2++) {
            auto partition1 = *cluster1->begin();
            auto partition2 = *cluster2->begin();
            ofs << i << " " << j << " " << wpJaccardDist(partition1, partition2) << endl;
            j++;
        }
        i++;
    }

    ofs.close();
    cout << "done!" << endl;
}

