#include <iostream>
#include <cstdio>
#include <vector>
#include <memory>
#include <cstring>

#include "joiner.h"

using namespace std;

Joiner::Joiner(int _numSeq, double** _score, const vector<string> & _proteinNames) {
  verbose = true;
  numSeq = _numSeq;
  score = new double*[numSeq];
  for (int i = 0; i < numSeq; i++) {
    score[i] = new double[numSeq];
    for (int j = 0; j < numSeq; j++) {
      score[i][j] = _score[i][j];
    }
  }
  proteinNames = _proteinNames;
  dist = new double*[numSeq];
  sumDist = new double*[numSeq];
  joinScore = new double*[numSeq];
  for (int i = 0; i < numSeq; i++) {
    dist[i] = new double[numSeq];
    sumDist[i] = new double[numSeq];
    joinScore[i] = new double[numSeq];
  }
  if (verbose) {
    printf("Number of sequences %d\n", numSeq);
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
      delete[] sumDist[i];
      delete[] joinScore[i];
    }
  delete[] score;
  delete[] dist;
  delete[] sumDist;
  delete[] joinScore;
}

void Joiner::join() {
  if (verbose) {
    printf("Start joining...\n");
  } 
  unique_ptr<bool[]> isMerged(new bool[numSeq]);
  for (int i = 0; i < numSeq; i++) {
    isMerged[i] = false;
  }
  for (int t = 0; t < numSeq - 1; t++) {
    updateDistance();
    if (verbose) {
      printf("\nDistance matrix \n");
      for (int i = 0; i < numSeq; i++) {
        for (int j = 0; j < numSeq; j++) {
          if (isMerged[i] || isMerged[j]) {
            continue;
          }
          printf("%.1lf ", dist[i][j]);
        }
        printf("\n");
      }
    }
    double bestScore = 1e9;
    int bestI = -1;
    int bestJ = -1;
    for (int i = 0; i < numSeq; i++) {
      if (!isMerged[i]) {
        for (int j = i + 1; j < numSeq; j++) {
          if (!isMerged[j]) {
            double curScore = (numSeq - 2) * dist[i][j] - sumDist[i][j];
            cout << i << ' ' << j << ' ' << curScore << endl;
            if (curScore < bestScore) {
              bestScore = curScore;
              bestI = i;
              bestJ = j;
            }
          }
        }
      }
    }
    if (verbose) {
      printf("BestScore = %.1lf\n", bestScore);
      printf("Merge %s into %s\n", proteinNames[bestJ].c_str(), proteinNames[bestI].c_str());
    }
    isMerged[bestJ] = true;
    parent[bestJ] = bestI;
    double sumI = 0;
    double sumJ = 0;
    for (int k = 0; k < numSeq; k++) {
      if (k == bestI || k == bestJ) {
        continue;
      }
      sumI += dist[bestI][k];
      sumJ += dist[bestJ][k];
    }
    for (int k = 0; k < numSeq; k++) {
      if (isMerged[k]) {
        continue;
      }
      if (k == bestI) {
        score[bestI][k] = (dist[bestI][bestJ] + 1.0 / (numSeq - 2) * (sumI - sumJ)) / 2;
      } else if (k == bestJ) {
        score[bestI][k] = score[k][bestI] = (dist[bestI][bestJ] + 1.0 / (numSeq - 2) * (sumJ - sumI)) / 2;
      } else {
        score[bestI][k] = score[k][bestI] = (dist[bestI][k] + dist[bestJ][k] - dist[bestI][bestJ]) / 2;
      }
    }
  }
}

void Joiner::updateDistance() {
  for (int i = 0; i < numSeq; i++) {
    for (int j = 0; j < numSeq; j++) {
      dist[i][j] = score[i][j];
      sumDist[i][j] = 0;
      for (int k = 0; k < numSeq; k++) {
        if (k != i && k != j) {
          sumDist[i][j] += score[i][k] + score[k][j];
        }
      }
    }
  }
}

void Joiner::setVerbose(bool value) {
  verbose = value;
}
