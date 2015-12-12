#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <memory>
#include <cstring>
#include <string>

#include "joiner.h"

using namespace std;

// Constructor
Joiner::Joiner(int _numSeq, double** _score, const vector<string> & _proteinNames) {
  verbose = false;
  numSeq = _numSeq;
  score = new double*[numSeq];
  for (int i = 0; i < numSeq; i++) {
    score[i] = new double[numSeq];
    for (int j = 0; j < numSeq; j++) {
      score[i][j] = _score[i][j];
    }
  }
  numNewNodes = 0;

  proteinNames = _proteinNames;
  dist = new double*[numSeq];
  sumDist = new double[numSeq];
  isMerged = new bool[numSeq];
  for (int i = 0; i < numSeq; i++) {
    dist[i] = new double[numSeq];
    isMerged[i] = false;
  }
  if (verbose) {
    printf("Number of sequences: %d\n", numSeq);
    printf("Protein names: ");
    for (int i = 0; i < numSeq; i++) {
      printf("%s ", proteinNames[i].c_str());
    }
    printf("\n");
    printf("Score matrix:\n");
    for (int i = 0; i < numSeq; i++) {
      for (int j = 0; j < numSeq; j++) {
        printf("%.1lf ", score[i][j]);
      }
      printf("\n");
    }
  }
}

Joiner::~Joiner() {
  for (int i = 0; i < numSeq; i++) {
      delete[] score[i];
      delete[] dist[i];
    }
  delete[] score;
  delete[] dist;
  delete[] sumDist;
}

// Neighbor-joining algorithm's main function
void Joiner::join() {
  if (verbose) {
    cout << endl;
    cout << "Start joining...\n" << endl;
  } 
  // Set all nodes to be unmerged
  for (int i = 0; i < numSeq; i++) {
    isMerged[i] = false;
  }
  int newNodeIndex = -1;
  // Loop through numSeq - 2 times
  for (int t = 0; t < numSeq - 2; t++) {
    // Reset the distance matrix 
    for (int i = 0; i < numSeq; i++) {
      for (int j = 0; j < numSeq; j++) {
        dist[i][j] = score[i][j];
      }
    }
    // Compute sum_k (dist[i][k]) for all i
    for (int i = 0; i < numSeq; i++) {
      sumDist[i] = 0;
      for (int k = 0; k < numSeq; k++) {
        if (!isMerged[k] && k != i) {
          sumDist[i] += dist[i][k];
        }
      }
    }
    if (verbose) {
      cout << "Distance matrix:" << endl;
      for (int i = 0; i < numSeq; i++) {
        if (isMerged[i]) {
          continue;
        }
        for (int j = 0; j < numSeq; j++) {
          if (isMerged[j]) {
            continue;
          }
          cout << "(" << i << "," << j << ")=" << dist[i][j] << "  ";
        }
        cout << "Sum = " << sumDist[i] << endl;
      }
    }
    double bestScore = 1e9;
    int bestI = -1;
    int bestJ = -1;
    if (verbose) {
      cout << "Join score matrix:" << endl;
    }
    // Compute the branch length matrix and find the maximum entry (bestI, bestJ)
    for (int i = 0; i < numSeq; i++) {
      if (!isMerged[i]) {
        for (int j = i + 1; j < numSeq; j++) {
          if (!isMerged[j]) {
            double curScore = (numSeq - 2 - t) * dist[i][j] - sumDist[i] - sumDist[j];
            if (curScore < bestScore) {
              bestScore = curScore;
              bestI = i;
              bestJ = j;
            }
            if (verbose) {
              cout << "(" << i << "," << j << ")=" << curScore << "  ";
            }
          }
        }
        if (verbose) {
          cout << endl;
        }
      }
    }
    // Make a new tree node
    string newNode = makeNewName();
    newNodeIndex = bestI;
    if (verbose) {
      cout << "Best score = " << bestScore << endl;
      cout << "Merge " << proteinNames[bestI] << " and " << proteinNames[bestJ] << " into " << newNode << endl;
      cout << endl;
    }
    // Re-compute the distance matrix
    for (int k = 0; k < numSeq; k++) {
      if (isMerged[k]) {
        continue;
      }
      if (k == bestI) {
        score[bestI][k] = (dist[bestI][bestJ] + 1.0 / (numSeq - 2 - t) * (sumDist[bestI] - sumDist[bestJ])) / 2;
      } else if (k == bestJ) {
        score[bestI][k] = score[k][bestI] = dist[bestI][bestJ] - score[bestI][bestI];
      } else {
        score[bestI][k] = score[k][bestI] = (dist[bestI][k] + dist[bestJ][k] - dist[bestI][bestJ]) / 2;
      }
    }
    // Add edges (bestI, newNode) and (bestJ, newNode) 
    addEdge(proteinNames[bestI], newNode, score[bestI][bestI]);
    addEdge(proteinNames[bestJ], newNode, score[bestI][bestJ]);
    proteinNames[bestI] = newNode;
    // Mark bestJ as having been merged
    isMerged[bestJ] = true;
    // Set distance from the new node to itself to be 0
    score[bestI][bestI] = 0;
  }
  // Merge the last sequence with last new node
  for (int i = 0; i < numSeq; i++) {
    if (!isMerged[i] && proteinNames[i] != proteinNames[newNodeIndex]) {
      addEdge(proteinNames[i], proteinNames[newNodeIndex], score[newNodeIndex][i]);
      break;
    }
  }
}

void Joiner::setVerbose(bool value) {
  verbose = value;
}

string Joiner::makeNewName() {
  return "New_" + to_string(numNewNodes++); 
}

void Joiner::addEdge(string x, string y, double score) {
  edges.push_back(make_pair(make_pair(x, y), score));
}

void Joiner::saveTree(ofstream &stream) {
  stream << edges.size() << endl;
  stream << "New_" << to_string(numNewNodes - 1) << endl;
  for (auto e : edges) {
    stream << e.first.first << ' ' << e.first.second << ' ' << e.second << endl;
  }
}
