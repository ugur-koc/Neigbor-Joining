#include <iostream>
#include <fstream>
#include <vector>

#include "joiner.h"

using namespace std;


void experiment1() {
  ifstream file;
  file.open("data2.txt");
  int numSeq;
  file >> numSeq;
  vector<string> proteinNames;
  for (int i = 0; i < numSeq; i++) {
    string s;
    file >> s;
    proteinNames.push_back(s);
  }
  double **score;
  score = new double*[numSeq];
  for (int i = 0; i < numSeq; i++) {
    score[i] = new double[numSeq];
    for (int j = 0; j < numSeq; j++) {
      file >> score[i][j];
    }
  }
  Joiner joiner(numSeq, score, proteinNames);
  joiner.join();
  file.close();
}

int main() {
  experiment1();
  return 0;
}
