#ifndef _JOINER_H_
#define _JOINER_H_

using namespace std;

class Joiner {
 public:
  Joiner(int _numSeq, double** _score, const vector<string> & _proteinNames);
  ~Joiner();
  void join();
  void setVerbose(bool value);
 private:
  bool verbose;
  // Number of protein sequences
  int numSeq;
  // Score matrix
  double **score;
  // Name of proteins
  vector<string> proteinNames;
  // Distance matrix
  double **dist;
  // Sum distance matrix sumDist[i] = sum_k(dist[i][k]) for k != i
  double *sumDist;
  // Join score matrix
  double **joinScore;
  // Parent of each node
  int* parent;
  // Has the node been merged?
  bool* isMerged;
};

#endif
