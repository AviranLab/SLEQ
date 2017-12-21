//
//  SHAPE-MAP.hpp
//  SHAPE-MAP
//
//  Created by Hua Li on 5/18/17.
//  Copyright © 2017 Hua Li. All rights reserved.
//

#ifndef SHAPE_SEQ_hpp
#define SHAPE_SEQ_hpp

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


int ConstructDMSEQ_old(int RNALen, vector<vector<bool>>Folds,double etaest, vector<double>gamaest,vector<vector<double>>&DSMatrix);


int ConstructDMSEQ_new(int RNALen, vector<vector<bool>>Folds,double etaest, vector<double>gamaest,vector<vector<double>>&DSMatrix);
template<typename T>
double mean(vector<T>in);




#endif /* SHAPE_SEQ_hpp */
