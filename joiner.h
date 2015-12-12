#ifndef _JOINER_H_
#define _JOINER_H_

#include <fstream>

using namespace std;

class Joiner {
 public:
  Joiner(int _numSeq, double** _score, const vector<string> & _proteinNames);
  ~Joiner();
  // Neighbor joining algorithm's main function
  void join();
  // Use to decide whether to print debugging strings or not
  void setVerbose(bool value);
  // Save phylogenetic tree to file
  void saveTree(ofstream &stream);
 private:
  bool verbose;
  // Number of protein sequences
  int numSeq;
  // Number of new nodes
  int numNewNodes;
  // Score matrix
  double **score;
  // Name of proteins
  vector<string> proteinNames;
  // Distance matrix
  double **dist;
  // Sum distance matrix sumDist[i] = sum_k(dist[i][k]) for k != i
  double *sumDist;
  // Has the node been merged?
  bool* isMerged;
  // Edges of the phylogenetic tree
  vector<pair<pair<string, string>, double>> edges;
  // Make new node name
  string makeNewName();
  // Add an edge (x, y) with weight "cost" to the phylogenetic tree
  void addEdge(string x, string y, double cost);
};

#endif
