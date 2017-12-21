//  Created by Hua Li
//  Copyright © 2016 Hua Li. All rights reserved.

#ifndef DMSTruncation_hpp
#define DMSTruncation_hpp

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


struct Pattern
{
    vector<bool>pat;
    int count;
    vector<int>muidx;
};



int ConstructFolds(string FileName, int RNALen, int SampleSize, int st, int ed, vector<vector<bool>>&Folds, int&FoldsNum);



int ConstructDMSEQ(int RNALen, vector<vector<bool>>FoldsA,vector<vector<bool>>FoldsC,double EtaA, double EtaC, vector<vector<double>>&DSMatrix);



int DrawReadsSEQ(int RNALen,int ReadsNum, vector<double>RandUnifData,vector<double>&ProbPlus, vector<int>& PlusReads);




template<typename T>
int runif(int len, T min, T max, vector<T>&out);

#endif /* DMSTruncation_hpp */
