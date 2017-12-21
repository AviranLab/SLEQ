//
//  Simulation.cpp
//  Simulation
//
//  Created by Hua Li on 5/19/17.
//  Copyright © 2017 Hua Li. All rights reserved.
//

#include "Simulation.hpp"

using namespace std;
vector<int> idx;
vector<int>pat;
int RNALen = 25;//change
int ReadsNum = 10000000;
int FoldsNum = 0;
vector<double>RelaAbund;
vector<double>RandUnifData(ReadsNum,0);
vector<int>MAPPlusReads;
vector<int>MAPMinusReads;
vector<double>PatternProbPlus;
vector<double>PatternProbMinus;
vector<vector<int>>Patterns;
vector<int>Pattern;
vector<vector<bool>>Folds;
double eta = 0.03;//change
vector<double>gama(RNALen,0);
string WritePath = "/Users/huali/Dropbox/";

double CurrentBinPlus = 0.0;
double CurrentBinMinus = 0.0;
int SearchBeginPlus = 0;
int SearchBeginMinus = 0;
bool stopflag = 0;
int t = 0;
double ThP = 0.99999; //change
double ThM = 0.99999;
int NoiseLevel = 1;//70dB//change

string RNAname = "BST";//change



int main(int argc, const char * argv[]) {
    ofstream myfilew;
    ifstream myfile;
    
    
    int SampleSize = 2;
    
    double rho[]={0.7,0.3};//change
    
    runif(RNALen,0.002,0.006,gama);
    
    vector<vector<bool>> DomiFolds = {//change
        //BST 0.5,0.5
          {0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0}
        //VcQrr3 0.15,0.4,0.2,0.25
        //{0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1},
        //{0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1},
        //{0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1},
        //{0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1}
        
        //TBWN 0.56,0.27,0.17
        /*
         {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0},
         {1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
         {1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,0,1,1,1,1,1,1,0,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,0,1,1}
         */
        //PreQ1 0.5,0.5  / 0.8,0.2
        
        //    {0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1},
        //  {0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1}
        
        //MST 0.5,0.3,0.2
        /*
         {1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,1,0,0,0,0,1,0,1,1,1,1,1,1,0,1,0,0,0,0,1},
         {1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,1,1},
         {1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,1}
         */
    };
    
    int DomiFoldsNum = (int)DomiFolds.size();
    
    //Read files & get unique folds & get indices of dominant folds
    vector<int>DomiFoldsIdx(DomiFoldsNum,0);
    string FilePath = "/Users/huali/Dropbox/New_Results_Review1/";
    string FileName = FilePath + RNAname+".txt";
    ConstructFolds(FileName, RNALen, SampleSize, DomiFolds, Folds, FoldsNum, DomiFoldsIdx);
    
    //Set relative abundances
    RelaAbund.assign(FoldsNum,0);
    
    //Add other minor folds
    vector<double>Otheridx(10,0);
    //
    /*
     runif(10,1.0,double(FoldsNum),Otheridx);
     vector<double>otherRela(10,0);
     runif(10,0.0,0.1,otherRela);
     for(int i = 0; i<10; ++i)
     RelaAbund[round(Otheridx[i])] = otherRela[i];
     for(int i = 0; i<DomiFoldsNum; ++i)
     RelaAbund[DomiFoldsIdx[i]] = 0;
     double sum = accumulate(RelaAbund.begin(), RelaAbund.end(), 0.0);
     for(int i = 0; i<FoldsNum; ++i)
     RelaAbund[i] = RelaAbund[i]/sum*0.1;
     for(int i = 0; i<10; ++i)
     cout<<round(Otheridx[i])<<": "<<RelaAbund[round(Otheridx[i])]<<endl;
     */
    //
    for(int i = 0; i<DomiFoldsNum; ++i)
        RelaAbund[DomiFoldsIdx[i]] = rho[i];
        
    double sum = accumulate(RelaAbund.begin(), RelaAbund.end(), 0.0);
    if(sum <0.9999 || sum>1.000001)
        cout<<"Error!";
    
    
    //Draw Reads for MAP
    runif(ReadsNum,0.0,1.0,RandUnifData);
    sort(RandUnifData.begin(),RandUnifData.end());
    for (int i = 0; i < RNALen; ++i) { idx.push_back(i); }
    DrawReadsMAPOrder(RNALen, ReadsNum, Folds, RelaAbund, eta, gama, RandUnifData, MAPPlusReads, MAPMinusReads, PatternProbPlus, PatternProbMinus, Patterns);
    int MAPPlusReadsSum = 0;
    int MAPMinusReadsSum = 0;
    double PatternProbPlusSum = 0.0;
    double PatternProbMinusSum = 0.0;
    MAPPlusReadsSum = accumulate(MAPPlusReads.begin(),MAPPlusReads.end(),0);
    MAPMinusReadsSum = accumulate(MAPMinusReads.begin(),MAPMinusReads.end(),0);
    PatternProbPlusSum = accumulate(PatternProbPlus.begin(),PatternProbPlus.end(),0.0);
    PatternProbMinusSum = accumulate(PatternProbMinus.begin(),PatternProbMinus.end(),0.0);
    int PatternsNum = (int)Patterns.size();
    MAPPlusReadsSum = accumulate(MAPPlusReads.begin(),MAPPlusReads.end(),0);
    MAPMinusReadsSum = accumulate(MAPMinusReads.begin(),MAPMinusReads.end(),0);
    for(int i = 0; i<PatternsNum; ++i)
    {
        MAPPlusReads[i] = round(double(MAPPlusReads[i])/double(MAPPlusReadsSum)*ReadsNum);
        MAPMinusReads[i] = round(double(MAPMinusReads[i])/double(MAPMinusReadsSum)*ReadsNum);
    }
    MAPPlusReadsSum = accumulate(MAPPlusReads.begin(),MAPPlusReads.end(),0);
    MAPMinusReadsSum = accumulate(MAPMinusReads.begin(),MAPMinusReads.end(),0);
    Trim(MAPPlusReads);
    Trim(MAPMinusReads);
    MAPPlusReadsSum = accumulate(MAPPlusReads.begin(),MAPPlusReads.end(),0);
    MAPMinusReadsSum = accumulate(MAPMinusReads.begin(),MAPMinusReads.end(),0);
    
    
    //Draw Reads for SEQ
    vector<int>SEQPlusReads(RNALen+1,0);
    vector<int>SEQMinusReads(RNALen+1,0);
    vector<double>SEQProbPlus(RNALen+1,0);
    vector<double>SEQProbMinus(RNALen+1,0);
    DrawReadsSEQ(RNALen,ReadsNum, Folds, RelaAbund, eta, gama, RandUnifData,SEQProbPlus,SEQProbMinus,SEQPlusReads,SEQMinusReads);
    int SEQPlusReadsSum = 0;
    int SEQMinusReadsSum = 0;
    double SEQProbPlusSum = 0.0;
    double SEQProbMinusSum = 0.0;
    SEQPlusReadsSum = accumulate(SEQPlusReads.begin(),SEQPlusReads.end(),0);
    SEQMinusReadsSum = accumulate(SEQMinusReads.begin(),SEQMinusReads.end(),0);
    SEQProbPlusSum = accumulate(SEQProbPlus.begin(),SEQProbPlus.end(),0.0);
    SEQProbMinusSum = accumulate(SEQProbMinus.begin(),SEQProbMinus.end(),0.0);
    
    //Add Noise
    //AddNoise(MAPPlusReads,ReadsNum);
    //MAPPlusReadsSum = accumulate(MAPPlusReads.begin(),MAPPlusReads.end(),0);
    //AddNoise(MAPMinusReads,ReadsNum);
    //MAPMinusReadsSum = accumulate(MAPMinusReads.begin(),MAPMinusReads.end(),0);
    //AddNoise(SEQPlusReads,ReadsNum);
    //SEQPlusReadsSum = accumulate(SEQPlusReads.begin(),SEQPlusReads.end(),0);
    //AddNoise(SEQMinusReads,ReadsNum);
    //SEQMinusReadsSum = accumulate(SEQMinusReads.begin(),SEQMinusReads.end(),0);
    
    //Construct MAP Profile (with noise)
    vector<double>MAPProfile(RNALen,0);
    vector<int>MAPHighReactIdx;
    vector<bool>MAPFoldsTemplate(RNALen,0);
    ConstructMAPProfile(RNALen, Patterns, MAPPlusReads, MAPMinusReads, MAPProfile, MAPHighReactIdx, MAPFoldsTemplate);
    
    

    //Kick out folds in MAP
    //vector<vector<bool>> MAPFoldsReduce;
    //vector<int>MAPFoldsReduceIdx;
    //int MAP_th = 0;
    //KickoutFolds(RNALen,MAP_th,Folds, MAPFoldsTemplate, MAPFoldsReduce,MAPFoldsReduceIdx);
    vector<vector<bool>> MAPFoldsReduce=Folds;

    
    //eta and gama estimation in MAP
    vector<double>MAPbetaest(RNALen,0);
    vector<double>MAPgamaest(RNALen,0);
    MAPParaEst(RNALen, ReadsNum, Patterns, MAPPlusReads, MAPMinusReads, MAPbetaest, MAPgamaest);
    for(int i = 0; i<RNALen; ++i)
    {
        cout<<gama[i]<<","<<MAPgamaest[i]<<endl;
    }
    
    //eta estimation in MAP
    double MAPetaest = 0.0;
    EtaEst(RNALen, MAPFoldsReduce, MAPbetaest, MAPHighReactIdx, MAPetaest);
    
    
    //beta and gama estimation in SEQ
    vector<double>SEQbetaest(RNALen,0);
    vector<double>SEQgamaest(RNALen,0);
    SEQParaEst(RNALen, SEQPlusReads, SEQMinusReads, SEQbetaest, SEQgamaest);
    
    //Construct SEQ Profile based on betaest (not reads due to signal decay);
    vector<double>SEQProfile(RNALen,0);
    vector<int>SEQHighReactIdx;
    vector<bool>SEQFoldsTemplate(RNALen,0);
    ConstructSEQProfile(RNALen, SEQbetaest,SEQProfile,SEQHighReactIdx,MAPFoldsTemplate);
    
    //Kick out folds in SEQ
    vector<vector<bool>>SEQFoldsReduce=Folds;
    //vector<int>SEQFoldsReduceIdx;
    //int SEQ_th = 0;
    //KickoutFolds(RNALen, SEQ_th,Folds, SEQFoldsTemplate, SEQFoldsReduce,SEQFoldsReduceIdx);
    
    //eta estimation in SEQ
    double SEQetaest = 0.0;
    EtaEst(RNALen, SEQFoldsReduce, SEQbetaest, SEQHighReactIdx, SEQetaest);
   
    //Coxnstruct Design Matrix in MAP
    int MAPFoldsReduceNum = (int)MAPFoldsReduce.size();
    vector<vector<double>>MAPDSMatrix(PatternsNum,vector<double>(MAPFoldsReduceNum,0));
    ConstructDMMAP(RNALen, MAPFoldsReduce, eta,MAPgamaest, Patterns,MAPDSMatrix);
    //ConstructDMMAP(RNALen, MAPFoldsReduce, eta,gama, Patterns,MAPDSMatrix);
    
    myfilew.open(WritePath+RNAname+"MAPDSMatrix.txt");
    for(int i = 0; i<PatternsNum; ++i)
        for(int j = 0; j<MAPFoldsReduceNum; ++j)
            myfilew<<MAPDSMatrix[i][j]<<"    ";
    myfilew.close();
    
    
    //Construct Design Matrix in SEQ
    int SEQFoldsReduceNum = (int)SEQFoldsReduce.size();
    vector<vector<double>>SEQDSMatrix(RNALen+1,vector<double>(SEQFoldsReduceNum,0));
    ConstructDMSEQ(RNALen, SEQFoldsReduce, SEQetaest, SEQgamaest, SEQDSMatrix);
    //ConstructDMSEQ(RNALen, SEQFoldsReduce, eta, gama, SEQDSMatrix);
    
    //Response in MAP and SEQ
    vector<double>MAPResponse(PatternsNum,0);
    vector<double>SEQResponse(RNALen+1,0);
    for(int i = 0; i<PatternsNum; ++i)
        MAPResponse[i] = double(MAPPlusReads[i])/ReadsNum;
    for(int i = 0; i<RNALen+1; ++i)
        SEQResponse[i] = double(SEQPlusReads[i])/ReadsNum;
    
    //    ofstream myfile;
    myfilew.open (WritePath+RNAname+"SEQResponse.txt");
    for(int i = 0; i<RNALen+1;++i)
        myfilew<<SEQResponse[i]<<"   ";
    myfilew.close();
    
    
    //Write SEQ DesignMatrix
    myfilew.open(WritePath+RNAname+"SEQDSMatrix.txt");
    for(int i = 0; i<RNALen+1; ++i)
        for(int j = 0; j<SEQFoldsReduceNum; ++j)
            myfilew<<SEQDSMatrix[i][j]<<"    ";
    myfilew.close();
    
    //Write MAP response
    myfilew.open (WritePath+RNAname+"MAPResponse.txt");
    for(int i = 0; i<PatternsNum;++i)
        myfilew<<MAPResponse[i]<<"   ";
    myfilew.close();
    
    //Write MAP DesignMatrix
    
    
    //Write SEQ files
    myfilew.open(WritePath+RNAname+"SEQpara.txt");
    myfilew<<RNALen<<"   ";
    myfilew<<SEQFoldsReduceNum<<"    ";
    for(int i = 0; i<SEQFoldsReduceNum; ++i)
        //myfilew<<SEQFoldsReduceIdx[i]+1<<" "; myfilew<<i+1<<" ";//without kicking
        myfilew<<i+1<<" ";//without kicking
    myfilew.close();
    
    //Write MAP files
    myfilew.open(WritePath+RNAname+"MAPpara.txt");
    myfilew<<RNALen<<"   ";
    myfilew<<MAPFoldsReduceNum<<"    ";
    myfilew<<PatternsNum<<"  ";
    for(int i = 0; i<MAPFoldsReduceNum; ++i)
        //myfilew<<MAPFoldsReduceIdx[i]+1<<" "; //with kicking
        myfilew<<i+1<<" ";//without kicking
    myfilew.close();
    return 0;
    
}

int ConstructDMSEQ(int RNALen, vector<vector<bool>>Folds,double etaest, vector<double>gamaest,vector<vector<double>>&DSMatrix)
{
    int FoldsNum = (int)Folds.size();
    for (int k = 0; k < RNALen+1; ++k)
    {
        for (int c = 0; c < FoldsNum; ++c)
        {
            float tempsum_plus = 1;
            if(k < RNALen) //Prob(adduct on site k)
            {
                if(Folds[c][k] == 1)
                    tempsum_plus = 1-(1-gamaest[k])*(1-etaest);
                else
                    tempsum_plus = gamaest[k];
            }
            for (int i = 0; i < k; ++i) //Prob(no adducts on sites from 0 to k-1)
            {
                if (Folds[c][i] == 1)
                    tempsum_plus = tempsum_plus*(1 - etaest)*(1 - gamaest[i]);
                else
                    tempsum_plus = tempsum_plus*(1 - gamaest[i]);
            }
            DSMatrix[k][c] = tempsum_plus;
        }
    }
    
    return 0;
}

int ConstructDMMAP(int RNALen, vector<vector<bool>>Folds, double etaest, vector<double>gamaest, vector<vector<int>>Patterns, vector<vector<double>>&DSMatrix)
{
    int FoldsNum = (int)Folds.size();
    int PatternsNum = (int)Patterns.size();
    for (int p = 0; p < PatternsNum; ++p)
    {
        vector<bool>pat(RNALen,0);
        for(vector<int>::iterator it = Patterns[p].begin(); it != Patterns[p].end(); ++it)
            pat[*it] = 1;
        for (int c = 0; c < FoldsNum; ++c)
        {
            float tempsum_plus = 1;
            for (int site = 0; site < RNALen; ++site)
            {
                if (pat[site] == 0)
                {
                    if (Folds[c][site] == 1)
                        tempsum_plus = tempsum_plus*(1 - etaest)*(1 - gamaest[site]);
                    else
                        tempsum_plus = tempsum_plus*(1 - gamaest[site]);
                }
                else
                {
                    if (Folds[c][site] == 1)
                        tempsum_plus = tempsum_plus*(1 - (1 - etaest)*(1 - gamaest[site]));
                    else
                        tempsum_plus = tempsum_plus*gamaest[site];
                }
            }
            DSMatrix[p][c] = tempsum_plus;
        }
    }
    return 0;
}

int EtaEst(int RNALen,vector<vector<bool>>Folds,vector<double>betaest,vector<int>HighReactIdx,double &etaest)
{
    vector<int> c;
    int FoldsNum = (int)Folds.size();
    for(int i = 0; i<RNALen; ++i)
    {
        for(int j = 0; j<FoldsNum; ++j)
        {
            if(Folds[j][i]==0)
                break;
            if(j == FoldsNum-1)
            {
                etaest += betaest[i];
                c.push_back(i);
            }
        }
        
    }
    vector<int>r;
    set_intersection(HighReactIdx.begin(), HighReactIdx.end(), c.begin(), c.end(), back_inserter(r));
    
    double s = 0.0;
    if(r.size()==0)
    {
        int len = (int)HighReactIdx.size();
        for(vector<int>::iterator it = HighReactIdx.begin(); it!=HighReactIdx.end(); ++it)
        {
            s += betaest[*it];
        }
        etaest = s/len;
    }
    else
    {
        int len = (int)r.size();
        for(vector<int>::iterator it = r.begin(); it!=r.end(); ++it)
        {
            s += betaest[*it];
        }
        etaest = s/len;
    }
    
    return 0;
}

int SEQParaEst(int RNALen, vector<int>SEQPlusReads, vector<int>SEQMinusReads, vector<double>&betaest, vector<double>&gamaest)
{
    for(int i = 0; i<RNALen; ++i)
    {
        double tempX = (double)accumulate(SEQPlusReads.begin()+i, SEQPlusReads.end(), 0);
        double tempY = (double)accumulate(SEQMinusReads.begin()+i, SEQMinusReads.end(), 0);
        gamaest[i] = double(SEQMinusReads[i])/tempY;
        betaest[i] = (double(SEQPlusReads[i])/tempX-gamaest[i])/(1-gamaest[i]);
        if(gamaest[i]<0) gamaest[i] = 0;
        if(betaest[i]<0) betaest[i] = 0;
    }
    
    return 0;
}

int MAPParaEst(int RNALen, int ReadsNum,vector<vector<int>>Patterns,vector<int>MAPPlusReads, vector<int>MAPMinusReads,vector<double> &betaest, vector<double>&gamaest)
{
    int PatternsNum = (int)Patterns.size();
    vector<int>mJc(RNALen,0);
    vector<int>pJc(RNALen,0);
    for(int p = 0; p<PatternsNum; ++p)
    {
        for(vector<int>::iterator it = Patterns[p].begin(); it!=Patterns[p].end(); ++it)
        {
            mJc[*it] += MAPMinusReads[p];
            pJc[*it] += MAPPlusReads[p];
        }
    }
    for(int i = 0; i<RNALen; ++i)
    {
        gamaest[i] = double(mJc[i])/ReadsNum;
        betaest[i] = (double(pJc[i])/ReadsNum - gamaest[i])/(1-gamaest[i]);
        if(gamaest[i]<0) gamaest[i] = 0;
        if(betaest[i]<0) betaest[i] = 0;
    }
    
    return 0;
}

int KickoutFolds(int RNALen,int th,vector<vector<bool>>Folds,vector<bool>MAPProfileTemplate,vector<vector<bool>>&FoldsReduce,vector<int>&FoldsReduceIdx)
{
    int FoldsNum = (int)Folds.size();
    for(int i = 0; i<FoldsNum; ++i)
    {
        int c = 0;
        for(int j = 0; j<RNALen; ++j)
        {
            if(MAPProfileTemplate[j]<Folds[i][j])
                ++c;
        }
        if(c<=th)
        {
            FoldsReduce.push_back(Folds[i]);
            FoldsReduceIdx.push_back(i);
        }
    }
    
    return 0;
}

int ConstructSEQProfile(int RNALen, vector<double>betaest,vector<double>&profile,vector<int>&HighReactIdx,vector<bool>&FoldsTemplate)
{
    double MaxVal = 0.0;
    for(int i = 0; i<RNALen; ++i)
    {
        if(betaest[i]>MaxVal)
            MaxVal = betaest[i];
    }
    for(int i = 0; i<RNALen; ++i)
    {
        profile[i] = betaest[i]/MaxVal;
        if(profile[i]>0.98) HighReactIdx.push_back(i);
        if(profile[i]>0.001) FoldsTemplate[i] = 1;
    }
    
    return 0;
}

int ConstructMAPchannelProfile(vector<vector<int>>Patterns,vector<int>MAPReads,vector<double>&profile)
{
    int PatternsNum = (int)Patterns.size();
    for(int p = 0; p<PatternsNum; ++p)
    {
        for(vector<int>::iterator it = Patterns[p].begin(); it!=Patterns[p].end(); ++it)
            profile[*it] += MAPReads[p];
        
    }
    return 0;
}

int ConstructMAPProfile(int RNALen,vector<vector<int>>Patterns,vector<int>MAPPlusReads,vector<int>MAPMinsuReads,vector<double>&profile,vector<int>&HighReactIdx, vector<bool>&FoldsTemplate)
{
    int PatternsNum = (int)Patterns.size();
    for(int p = 0; p<PatternsNum; ++p)
    {
        for(vector<int>::iterator it = Patterns[p].begin(); it!=Patterns[p].end(); ++it)
            profile[*it] += (MAPPlusReads[p]-MAPMinsuReads[p]);
        /*        vector<bool>pat(RNALen,0);
         Double2Binary(Patterns[p], RNALen, pat);
         for(int site = 0; site<RNALen; ++site)
         {
         if(pat[site]==1)
         {
         profile[site] += (MAPPlusReads[p]-MAPMinsuReads[p]);
         }
         }
         */
    }
    double MaxVal = 0.0;
    for(int i = 0; i<RNALen; ++i)
    {
        if(profile[i]>MaxVal)
            MaxVal = profile[i];
    }
    for(int i = 0; i<RNALen; ++i)
    {
        profile[i] = profile[i]/MaxVal;
        if(profile[i]>=0.95)
            HighReactIdx.push_back(i);
        if(profile[i]>0.001)
            FoldsTemplate[i] = 1;
    }
    
    return 0;
}

int AddNoise(vector<int>&Signal,int ReadsNum)
{
    int len = (int)Signal.size();
    for(int i = 0; i<len; ++i)
    {
        double sigma_noise = Signal[i]/pow(10,NoiseLevel);
        unsigned seed = (unsigned)chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator(seed);
        normal_distribution<double> distribution(0.0,sigma_noise);
        int noise = (int)round(distribution(generator));
        Signal[i] += noise;
        if(Signal[i]<0)
            Signal[i] = 0;
    }
    int rem = accumulate(Signal.begin(), Signal.end(), 0)-ReadsNum;
    while(rem>0)
    {
        for(int i = 0; i<len&&rem>0; ++i)
        {
            if(Signal[i]>0)
            {
                Signal[i] -= 1;
                --rem;
            }
        }
    }
    while(rem<0)
    {
        for(int i = 0; i<len&&rem<0; ++i,++rem)
        {
            Signal[i] += 1;
        }
    }
    
    return 0;
}

int Trim(vector<int>&data)
{
    int rem = accumulate(data.begin(), data.end(), 0)-ReadsNum;
    int len = (int)data.size();
    while(rem>0)
    {
        for(int i = 0; i<len&&rem>0; ++i)
        {
            if(data[i]>0)
            {
                data[i] -= 1;
                --rem;
            }
        }
    }
    while(rem<0)
    {
        for(int i = 0; i<len&&rem<0; ++i,++rem)
        {
            data[i] += 1;
        }
    }
    return 0;
}

int DrawReadsSEQ(int RNALen,int ReadsNum, vector<vector<bool>>Folds, vector<double>RelaAbund, double eta, vector<double>gama, vector<double>RandUnifData,vector<double>&ProbPlus, vector<double>&ProbMinus, vector<int>&PlusReads,vector<int>&MinusReads)
{
    int FoldsNum = (int)Folds.size();
    for (int k = 0; k < RNALen+1; ++k)
    {
        ProbPlus[k] = 0;
        ProbMinus[k] = 0;
        for (int c = 0; c < FoldsNum; ++c)
        {
            if(RelaAbund[c]==0)
                continue;
            float tempsum_plus = 1;
            float tempsum_minus = 1;
            if(k < RNALen) //Prob(adduct on site k)
            {
                tempsum_minus = gama[k];
                if(Folds[c][k] == 1)
                    tempsum_plus = 1-(1-gama[k])*(1-eta);
                else
                    tempsum_plus = gama[k];
            }
            for (int i = 0; i < k; ++i) //Prob(no adducts on sites from 0 to k-1)
            {
                tempsum_minus = tempsum_minus*(1 - gama[i]);
                if (Folds[c][i] == 1)
                    tempsum_plus = tempsum_plus*(1 - eta)*(1 - gama[i]);
                else
                    tempsum_plus = tempsum_plus*(1 - gama[i]);
            }
            ProbPlus[k] += tempsum_plus*RelaAbund[c];
            ProbMinus[k] += tempsum_minus*RelaAbund[c];
        }
    }
    vector<double>ProbPlusCumsum(RNALen+1,0);
    vector<double>ProbMinusCumsum(RNALen+1,0);
    partial_sum(ProbPlus.begin(), ProbPlus.end(), ProbPlusCumsum.begin());
    partial_sum(ProbMinus.begin(), ProbMinus.end(), ProbMinusCumsum.begin());
    int PlusBegin = 0;
    int MinusBegin = 0;
    for(int i = 0; i<RNALen+1; ++i)
    {
        for(int j = PlusBegin; j<ReadsNum; ++j)
        {
            if(RandUnifData[j]>ProbPlusCumsum[i])
            {
                PlusReads[i] = j - PlusBegin;
                PlusBegin = j;
                break;
            }
            if(j == ReadsNum-1)
            {
                PlusReads[i] = ReadsNum - PlusBegin;
                PlusBegin = ReadsNum;
            }
        }
    }
    for(int i = 0; i<RNALen+1; ++i)
    {
        for(int j = MinusBegin; j<ReadsNum; ++j)
        {
            if(RandUnifData[j]>ProbMinusCumsum[i])
            {
                MinusReads[i] = j - MinusBegin;
                MinusBegin = j;
                break;
            }
            if(j == ReadsNum-1)
            {
                MinusReads[i] = ReadsNum - MinusBegin;
                MinusBegin = ReadsNum;
            }
        }
    }
    
    return 0;
}


int DrawReadsMAPOrder(int RNALen, int ReadsNum, vector<vector<bool>>Folds,vector<double>RelaAbund, double eta, vector<double>gama, vector<double>RandUnifData,vector<int>&PlusReads,vector<int>&MinusReads,vector<double>&PatternProbPlus,vector<double>&PatternProbMinus,vector<vector<int>>&Patterns)
{
    int i = 0;
    for(i = 0; i<RNALen && stopflag ==0; ++i)
    {
        go(0,i);
    }
    return 0;
}

int ConstructFolds(string FileName, int RNALen, int SampleSize, vector<vector<bool>> DomiFolds,vector<vector<bool>>&Folds, int&FoldsNum, vector<int>&DomiFoldsIdx)
{
    ifstream myfile;
    myfile.open (FileName.c_str());
    char * dump = new char[6];
    char * Buffer = new char[RNALen];
    vector<bool> FoldsUnit(RNALen,0);
    FoldsNum = 0;
    
    int dumpnum = 0;
    //Read Files, convert to 0,1 form, remain unique folds
    for(int i = 0; i<10; ++i)
    {
        myfile.read(dump,sizeof(char)*1);
        if(dump[0]=='('||dump[0]==')'||dump[0]=='.')
        {
            Buffer[0] = dump[0];
            myfile.read(&Buffer[1],sizeof(char)*(RNALen-1));
            if(i == 5) dumpnum = 1;
            if(i == 6) dumpnum = 2;
            break;
        }
    }
    myfile.read(dump,sizeof(char)*dumpnum);
    for(int j = 0; j<RNALen; ++j)
    {
        if(Buffer[j] == '.')
            FoldsUnit[j] = 1;
        else
            FoldsUnit[j] = 0;
    }
    Folds.push_back(FoldsUnit);
    ++FoldsNum;
    for(int i = 1; i<SampleSize; ++i)
    {
        bool flag = 0;
        myfile.read(Buffer,sizeof(char)*RNALen);
        myfile.read(dump,sizeof(char)*dumpnum);
        for(int j = 0; j<RNALen; ++j)
        {
            if(Buffer[j] == '.')
                FoldsUnit[j] = 1;
            else
                FoldsUnit[j] = 0;
        }
        for(int k = 0; k<FoldsNum; ++k)
        {
            if(Folds[k] == FoldsUnit)
            {
                
                flag = 1;
                break;
            }
        }
        if(flag == 0)
        {
            Folds.push_back(FoldsUnit);
            ++FoldsNum;
        }
    }
    myfile.close();
    delete []dump;
    delete []Buffer;
    //find indeces of dominant folds
    int DomiFoldsNum = (int)DomiFolds.size();
    for(int i = 0; i<DomiFoldsNum; ++i)
    {
        for(int j = 0; j<FoldsNum; ++j)
        {
            if(Folds[j] == DomiFolds[i])
            {
                DomiFoldsIdx[i] = j;
                break;
            }
        }
    }
    
    return 0;
}

int Int2Binary(int num,int len,vector<bool>&out)
{
    
    int count = len-1;
    while(num!=0 && count>=0)
    {
        out[count] = num%2;
        num /= 2;
        --count;
    }
    return 0;
}

int Double2Binary(double num,int len,vector<bool>&out)
{
    
    int count = len-1;
    while( num>=1 && count>=0)
    {
        if(num/2==round(num/2))
            out[count] = 0;
        else
            out[count] = 1;
        num = floor(num/2);
        count--;
    }
    return 0;
}

int Binary2BDouble(vector<bool>in,double &out)
{
    int len = (int)in.size();
    out = 0.0;
    for(int i = 0; i<len; ++i)
    {
        if(in[i]==1)
            out += pow(2,RNALen-1-i);
    }
    return 0;
}

template<typename T>
int rnorm(int len, T m, T sd,vector<T>&out)
{
    unsigned seed = (unsigned)chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<T> distribution(m,sd);
    
    for (int i=0; i<len; ++i) {
        out[i]= distribution(generator);
    }
    return 0;
}


template<typename T>
int runif(int len, T min, T max, vector<T>&out)
{
    unsigned seed = (unsigned)chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_real_distribution<T> distribution(min,max);
    for (int i=0; i<len; ++i) {
        out[i]= distribution(generator);
    }
    return 0;
}

template<typename T>
double mean(vector<T>in)
{
    int len = (int)in.size();
    double sum = accumulate(in.begin(), in.end(), 0.0);
    return sum/len;
}

template<typename T>
double var(vector<T>in)
{
    int len = (int)in.size();
    double m = mean(in);
    double s = 0;
    for(int i = 0; i<len; ++i)
        s+=(in[i]-m)*(in[i]-m);
    return s/(len-1);
}

void go(int offset, int k) {
    if(stopflag == 1)
        return;
    if (k == 0) {
        vector<bool>patunit(RNALen,0);
        for(vector<int>::iterator it = pat.begin(); it!=pat.end(); ++it)
            patunit[*it] = 1;
        double PatternProbPlusCurrent = 0.0;
        double PatternProbMinusCurrent = 0.0;
        for (int c = 0; c < FoldsNum; ++c)
        {
            if(RelaAbund[c]==0)
                continue;
            double tempsum_plus = 1.0;
            double tempsum_minus = 1.0;
            for (int site = 0; site < RNALen; ++site)
            {
                if (patunit[site] == 0)
                {
                    tempsum_minus = tempsum_minus*(1 - gama[site]);
                    if (Folds[c][site] == 1)
                        tempsum_plus = tempsum_plus*(1 - eta)*(1 - gama[site]);
                    else
                        tempsum_plus = tempsum_plus*(1 - gama[site]);
                }
                else
                {
                    tempsum_minus = tempsum_minus*gama[site];
                    if (Folds[c][site] == 1)
                        tempsum_plus = tempsum_plus*(1 - (1 - eta)*(1 - gama[site]));
                    else
                        tempsum_plus = tempsum_plus*gama[site];
                }
            }
            PatternProbPlusCurrent += tempsum_plus*RelaAbund[c];
            PatternProbMinusCurrent += tempsum_minus*RelaAbund[c];
        }
        
        int PlusReadsCurrent = 0;
        CurrentBinPlus += PatternProbPlusCurrent;
        
        for(int i = SearchBeginPlus; i<ReadsNum; ++i)
        {
            if(RandUnifData[i]>CurrentBinPlus)
            {
                PlusReadsCurrent=i-SearchBeginPlus;
                SearchBeginPlus = i;
                break;
            }
            
            if(i == ReadsNum-1)
            {
                PlusReadsCurrent = ReadsNum - SearchBeginPlus;
                SearchBeginPlus = ReadsNum;
            }
        }
        
        
        int MinusReadsCurrent = 0;
        CurrentBinMinus += PatternProbMinusCurrent;
        for(int i = SearchBeginMinus; i<ReadsNum; ++i)
        {
            if(RandUnifData[i]>CurrentBinMinus)
            {
                MinusReadsCurrent=i-SearchBeginMinus;
                SearchBeginMinus = i;
                break;
            }
            if(i == ReadsNum-1)
            {
                MinusReadsCurrent = ReadsNum - SearchBeginMinus;
                SearchBeginMinus = ReadsNum;
            }
        }
        
        if(PlusReadsCurrent!=0 || MinusReadsCurrent!=0)
        {
            //            double patnum = 0.0;
            //           Binary2BDouble(patunit,patnum);
            MAPPlusReads.push_back(PlusReadsCurrent);
            MAPMinusReads.push_back(MinusReadsCurrent);
            PatternProbPlus.push_back(PatternProbPlusCurrent);
            PatternProbMinus.push_back(PatternProbMinusCurrent);
            Patterns.push_back(pat);
            ++t;
        }
        if((SearchBeginPlus==ReadsNum&&SearchBeginMinus==ReadsNum)|| (CurrentBinPlus>ThP&&CurrentBinMinus>ThM))
        {
            stopflag = 1;
            return;
        }
        //        pretty_print(pattemp);
        return;
    }
    for (int i = offset; i <= idx.size() - k; ++i) {
        if(stopflag == 1)
            return;
        pat.push_back(idx[i]);
        go(i+1, k-1);
        pat.pop_back();
    }
}

