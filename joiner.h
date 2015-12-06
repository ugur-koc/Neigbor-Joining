#ifndef _JOINER_H_
#define _JOINER_H_

#include <fstream>

using namespace std;

class Joiner {
 public:
  Joiner(int _numSeq, double** _score, const vector<string> & _proteinNames);
  ~Joiner();
  void join();
  void setVerbose(bool value);
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
  // Make new node name
  vector<pair<pair<string, string>, double>> edges;
  string makeNewName();
  void addEdge(string x, string y, double cost);
};

#endif
