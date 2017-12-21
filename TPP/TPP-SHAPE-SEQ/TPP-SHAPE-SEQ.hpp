//
//  TPP-SHAPE-SEQ.hpp
//  ADD-SHAPE-SEQ
//
//  Created by Hua Li on 5/19/17.
//  Copyright © 2017 Hua Li. All rights reserved.
//

#ifndef TPP_SHAPE_SEQ_hpp
#define TPP_SHAPE_SEQ_hpp

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

int KickoutFolds(vector<int> openidx,vector<vector<bool>>Folds,vector<vector<bool>>&FoldsReduce,vector<int>&FoldsReduceIdx);



int ConstructFolds(string FileName, int RNALen, int SampleSize, vector<vector<bool>>&Folds, int&FoldsNum);


int ConstructDMSEQ(int RNALen, vector<vector<bool>>Folds,double etaest, vector<double>gamaest,vector<vector<double>>&DSMatrix);

template<typename T>
double mean(vector<T>in);


#endif /* TPP_SHAPE_SEQ_hpp */
