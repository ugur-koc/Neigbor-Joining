#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

#include "joiner.h"

using namespace std;

int numSeq;
vector<string> proteinNames;
double **score;

void loadData(ifstream &stream) {
  stream >> numSeq;
  for (int i = 0; i < numSeq; i++) {
    string tmp;
    stream >> tmp;
    proteinNames.push_back(tmp);
  }
  score = new double*[numSeq];
  for (int i = 0; i < numSeq; i++) {
    score[i] = new double[numSeq];
    for (int j = 0; j < numSeq; j++) {
      stream >> score[i][j];
    }
  }
}

void experiment1() {
  ifstream file("data2.txt");
  loadData(file);
  Joiner joiner(numSeq, score, proteinNames);
  joiner.join();
  file.close();
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    cout << endl;
    cout << "USAGE: ./main <INPUT> <OUTPUT> <VERBOSE>"  << endl;
    cout << "    Input: file containing protein names and a distance matrix" << endl;
    cout << "    Output: file to store the computed phylogenetic tree" << endl;
    cout << "    Verbose: \"verbose\" to print debugging logs" << endl;
    cout << endl;
    return 1;
  }
  string input = argv[1];
  string output = argv[2];
  
  ifstream infile(input);
  if (infile) {
    loadData(infile);
    infile.close();
  } else {
    cout << "Input not found!" << endl;
    return 1;
  }
  Joiner joiner(numSeq, score, proteinNames);
  if (argv[3] != NULL && strcmp(argv[3], "verbose") == 0) {
    joiner.setVerbose(true);
  }
  joiner.join();
  ofstream outfile(output);
  joiner.saveTree(outfile);
  outfile.close();
  return 0;
}
