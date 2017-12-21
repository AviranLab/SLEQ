//
//  Simulation.hpp
//  Simulation
//
//  Created by Hua Li on 5/19/17.
//  Copyright © 2017 Hua Li. All rights reserved.
//

#ifndef Simulation_hpp
#define Simulation_hpp

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

int EtaEst(int RNALen,vector<vector<bool>>Folds,vector<double>betaest,vector<int>HighReactIdx,double &etaest);

int SEQParaEst(int RNALen, vector<int>SEQPlusReads, vector<int>SEQMinusReads, vector<double>&betaest, vector<double>&gamaest);

int MAPParaEst(int RNALen, int ReadsNum,vector<vector<int>>Patterns,vector<int>MAPPlusReads, vector<int>MAPMinusReads,vector<double> &betaest, vector<double>&gamaest);

int KickoutFolds(int RNALen,int th,vector<vector<bool>>Folds,vector<bool>MAPProfileTemplate,vector<vector<bool>>&FoldsReduce,vector<int>&FoldsReduceIdx);

int ConstructSEQProfile(int RNALen, vector<double>betaest,vector<double>&profile,vector<int>&HighReactIdx,vector<bool>&FoldsTemplate);


int ConstructMAPchannelProfile(vector<vector<int>>Patterns,vector<int>MAPReads,vector<double>&profile);

int ConstructMAPProfile(int RNALen,vector<vector<int>>Patterns,vector<int>MAPPlusReads,vector<int>MAPMinsuReads,vector<double>&profile,vector<int>&HighReactIdx, vector<bool>&FoldsTemplate);

int AddNoise(vector<int>&Signal,int ReadsNum);

int Trim(vector<int>&data);

int ConstructFolds(string FileName, int RNALen, int SampleSize, vector<vector<bool>> DomiFolds,vector<vector<bool>>&Folds, int&FoldsNum, vector<int>&DomiFoldsIdx);

int DrawReadsMAP(int RNALen, int ReadsNum, vector<vector<bool>>Folds,vector<double>RelaAbund, double eta, vector<double>gama, vector<double>RandUnifData,vector<int>&PlusReads,vector<int>&MinusReads,vector<double>&PatternProbPlus,vector<double>&PatternProbMinus,vector<vector<int>>&Patterns);

int DrawReadsSEQ(int RNALen,int ReadsNum, vector<vector<bool>>Folds, vector<double>RelaAbund, double eta, vector<double>gama, vector<double>RandUnifData,vector<double>&ProbPlus, vector<double>&ProbMinus, vector<int>&PlusReads,vector<int>&MinusReads);

int ConstructDMMAP(int RNALen, vector<vector<bool>>Folds, double etaest, vector<double>gamaest, vector<vector<int>>Patterns, vector<vector<double>>&DSMatrix);

int ConstructDMSEQ(int RNALen, vector<vector<bool>>Folds,double etaest, vector<double>gamaest,vector<vector<double>>&DSMatrix);

int DrawReadsMAPOrder(int RNALen, int ReadsNum, vector<vector<bool>>Folds,vector<double>RelaAbund, double eta, vector<double>gama, vector<double>RandUnifData,vector<int>&PlusReads,vector<int>&MinusReads,vector<double>&PatternProbPlus,vector<double>&PatternProbMinus,vector<vector<int>>&Patterns);

//Small Functions
template<typename T>
int runif(int len, T min, T max, vector<T>&out);

template<typename T>
int rnorm(int len, T m, T sd,vector<T>&out);

int Int2Binary(int num,int len,vector<bool>&out);

int Double2Binary(double num,int len,vector<bool>&out);

int Binary2BDouble(vector<bool>in,double &out);

template<typename T>
double mean(vector<T>in);

template<typename T>
double var(vector<T>in);

void pretty_print(const vector<int>& v);

void go(int offset, int k);

#endif /* Simulation_hpp */
