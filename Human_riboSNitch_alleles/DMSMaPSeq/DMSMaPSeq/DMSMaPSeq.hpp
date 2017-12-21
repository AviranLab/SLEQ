//
//  RNAEnsembleRealMAP.hpp
//  RNAEnsembleRealMAP
//
//  Created by Hua Li on 12/10/16.
//  Copyright © 2016 Hua Li. All rights reserved.
//

#ifndef DMSMaPSeq_hpp
#define DMSMaPSeq_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include <valarray>
#include <numeric>
#include <ctime>
#include <chrono>


using namespace std;

struct sam
{
    char ref;
    int pos;
    vector<char> cigar;
    vector<char> seq;
};

struct Pattern
{
    vector<bool>pat;
    int count;
    vector<int>muidx;
};



int ConstructFolds(string FileName, int RNALen, int SampleSize, int st, int ed, vector<vector<bool>>&Folds, int&FoldsNum);


int ConstructDMMAP(int RNALen, vector<vector<bool>>FoldsA, vector<vector<bool>>FoldsC,vector<Pattern>ReadsAC, double EtaA, double EtaC,vector<vector<float>>&DSMatrix);

int DrawReadsSEQ(int RNALen,int ReadsNum, vector<double>RandUnifData,vector<double>&ProbPlus, vector<int>& PlusReads);
int DrawReadsMAP(int PatternsNum,int ReadsNum, vector<double>RandUnifData,vector<double>&ProbPlus, vector<int>& PlusReads);



template<typename T>
int runif(int len, T min, T max, vector<T>&out);

#endif /* DMSMaPSeq_hpp */
