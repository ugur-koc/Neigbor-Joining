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
  void updateDistance();
  bool verbose;
  // Number of protein sequences
  int numSeq;
  // Score matrix
  double** score;
  // Name of proteins
  vector<string> proteinNames;
  // Distance matrix
  double **dist;
  // Sum distance matrix sumDist[i][j] = sum_k(dist[i][k] + dist[k][j]) for k != i, j
  double **sumDist;
  // Join score matrix
  double **joinScore;
  // Parent of each node 
  int* parent;
};

#endif
