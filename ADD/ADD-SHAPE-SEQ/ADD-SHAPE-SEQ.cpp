//
//  ADD-SHAPE-SEQ.cpp
//  ADD-SHAPE-SEQ
//
//  Created by Hua Li on 5/18/17.
//  Copyright © 2017 Hua Li. All rights reserved.
//

#include "ADD-SHAPE-SEQ.hpp"

using namespace std;
string Folder = "/Users/huali/Desktop/Codes/ADD/"; //Set the path where the ensemble file locates. For example, ADD74.txt
string RNAname = "ADD74";



//ADD
int St = 1;//from 3'
int Ed = 70;//from 3' remove last 3 nucleotides in 5', as paper only gives 71 nt structures, so complete reads
int RNALenFold = 74;
int RNALenNNLS = Ed-St +1; //no complete reads


int main(int argc, const char * argv[]) {
    ofstream myfilew;
    ifstream myfile;
    
    //Read files & get unique folds & get indices of dominant folds
    int SampleSize = 1000;
    int FoldsNum = 0;
    vector<vector<bool>>Folds;//from 3' to 5'
    string FileName = Folder + RNAname+".txt";
    ConstructFolds(FileName, RNALenFold, SampleSize, Folds, FoldsNum);//Folds from 3' to 5'
    
    
    //Order: complete seq, 5' to 3'
    vector<int>SEQPlusReads=//change
    //ADD REP1 v2.1
     {12185,  3230, 47464,  2786,  1465,   416,   345,  4469,   330, 10911,
     2366,  3500,  9202, 12070,  7025,  3554,   177,   170,  1730,   800,
     2727,   501,   709,  6930,  2878,  2023,  1078,   528,  3492,   502,
     278,   127,   162,    82,   162,   430,  2438,   744,  2183, 23945,
     1113,  1606,  9074,  1241,   312,   535,   146,   112,    40,    95,
     89,   174,  1070,  4610,  1161,  1149,   302,    67,   244,   134,
     145,    69,    50,  1388,   781,   774,  2065,   112,    13,    71,
     73,    18,    65,    90};
    vector<int>SEQMinusReads=
    //ADD rep1 v2.1
     {8592, 1461, 8502,  253,  731,  355,  250, 4280,  265, 6719,
     1926, 2416, 4335, 1137,  583, 1194,  110,  147, 1392,  429,
     1682,  195,  217, 2529, 2172, 1506,  143,  287,  791,  183,
     82,   59,  102,   57,  129,  316,  848,  414,  265, 1656,
     788,  252, 2574,  468,  120,  265,   81,   51,   16,   52,
     37,   64,  299,  510,   88,  488,  121,   28,   77,   41,
     38,   22,   28,  455,   84,   91,  799,   52,    4,   26,
     27,    6,   74,  262};
    
    
    //Order: beta and gama here are from 3' to 5'
    vector<double>SEQbetaest(RNALenFold-1,0);//remove the last nucleotide
    vector<double>SEQgamaest(RNALenFold-1,0);
    vector<double>SEQRate_Plus(RNALenFold-1,0);
    for(int i = 0; i<RNALenFold-1; ++i)
    {
        //Must have complete reads
        
        double tempX = (double)accumulate(SEQPlusReads.begin(), SEQPlusReads.begin()+i+2, 0);
        double tempY = (double)accumulate(SEQMinusReads.begin(), SEQMinusReads.begin()+i+2, 0);
        SEQgamaest[RNALenFold-i-2] = double(SEQMinusReads[i+1])/tempY;
        SEQRate_Plus[RNALenFold-i-2] = double(SEQPlusReads[i+1])/tempX;
        SEQbetaest[RNALenFold-i-2] = (SEQRate_Plus[RNALenFold-i-2]-SEQgamaest[RNALenFold-i-2])/(1-SEQgamaest[RNALenFold-i-2]);
        if(SEQgamaest[RNALenFold-i-2]<0) SEQgamaest[RNALenFold-i-2]=0;
        if(SEQbetaest[RNALenFold-i-2]<0) SEQbetaest[RNALenFold-i-2]=0;
        
    }
    for(int i=Ed; i<RNALenFold-1; ++i) // set as 0, otherwise, it affects SEQHIghReactIdx, it may be in the first St sites.
    {
        SEQgamaest[i]=0.0;
        SEQbetaest[i]=0.0;
    }
    //ETA
    double SEQetaest = 0.0;
    double betamean = mean(SEQbetaest);
    vector<double>SEQbetasort(SEQbetaest);
    sort(SEQbetasort.begin(),SEQbetasort.end(),greater<double>());
    int stopidx = 0;
    for (int i = 0; i<RNALenFold-1; ++i)
    {
        if(SEQbetasort[i]<betamean)
        {
            stopidx = i;
            break;
        }
    }
    if(stopidx%2==0)
        SEQetaest = (SEQbetasort[stopidx*0.5]+SEQbetasort[stopidx*0.5-1])/2;
    else
        SEQetaest = SEQbetasort[(stopidx-1)*0.5];
    
  
    
    //Kick out folds in SEQ
    vector<vector<bool>>SEQFoldsReduce;
    vector<int>SEQFoldsReduceIdx;
    vector<int>opensite;
    vector<double>highbeta;
    for(int i = 0; i<RNALenFold-1; ++i)
        if(SEQbetaest[i]>SEQetaest)
        {
            opensite.push_back(i);
            highbeta.push_back(SEQbetaest[i]);
        }
    KickoutFolds(opensite,Folds, SEQFoldsReduce,SEQFoldsReduceIdx);//Ignore the last one, only compare the RNALenFold-1 to determine kickout.
    
    for(int i=Ed; i<RNALenFold-1; ++i)
    {
        SEQgamaest.pop_back();
        SEQbetaest.pop_back();
    }
    
    

    //Construct Design Matrix in SEQ
    int SEQFoldsReduceNum = (int)SEQFoldsReduce.size();
    vector<vector<double>>SEQDSMatrix(RNALenNNLS,vector<double>(SEQFoldsReduceNum,0));
    ConstructDMSEQ(Ed-St+1, SEQFoldsReduce, SEQetaest, SEQgamaest, SEQDSMatrix);
    
    
    //Response in SEQ
    vector<double>SEQResponse(RNALenNNLS,0);
    //from 3' to 5'
    for(int i = St; i<=Ed; ++i)
        SEQResponse[i-St] = SEQPlusReads[RNALenFold-i];//because plus reads already remove the last nt, so it does not need -1
    double ReadsNum = accumulate(SEQResponse.begin(), SEQResponse.end(), 0.0);
    for(int i = 0; i<RNALenNNLS; ++i)
        SEQResponse[i] = double(SEQResponse[i])/ReadsNum;
    

    //Write SEQ responses
    myfilew.open (Folder+RNAname+"SEQResponse.txt");
    for(int i = 0; i<RNALenNNLS;++i)
        myfilew<<SEQResponse[i]<<"   ";
    myfilew.close();
    
    
    //Write SEQ DesignMatrix
    myfilew.open(Folder+RNAname+"SEQDSMatrix.txt");
    for(int i = 0; i<RNALenNNLS; ++i)
        for(int j = 0; j<SEQFoldsReduceNum; ++j)
            myfilew<<SEQDSMatrix[i][j]<<"    ";
    myfilew.close();
    
    
    //Write SEQ files
    myfilew.open(Folder+RNAname+"SEQpara.txt");
    myfilew<<Ed-St+1<<"   ";
    myfilew<<SEQFoldsReduceNum<<"    ";
    for(int i = 0; i<SEQFoldsReduceNum; ++i)
        myfilew<<SEQFoldsReduceIdx[i]+1<<" ";
    myfilew.close();
    
    myfilew.open (Folder+RNAname+"SEQProfile.txt");
    for(int i = 0; i<RNALenNNLS;++i)
        myfilew<<SEQbetaest[i]<<"   ";
    myfilew.close();

    
    return 0;
}

int ConstructDMSEQ(int RNALen, vector<vector<bool>>Folds,double etaest, vector<double>gamaest,vector<vector<double>>&DSMatrix)
{
    int FoldsNum = (int)Folds.size();
    for (int k = 0; k < RNALen; ++k)
    {
        for (int c = 0; c < FoldsNum; ++c)
        {
            float tempsum_plus = 1;
            if(k < RNALen) //Prob(adduct on site k)
            {
                if(Folds[c][k+St] == 1)
                    tempsum_plus = 1-(1-gamaest[k])*(1-etaest);
                else
                    tempsum_plus = gamaest[k];
            }
            for (int i = 0; i < k; ++i) //Prob(no adducts on sites from 0 to k-1)
            {
                if (Folds[c][i+St] == 1)
                    tempsum_plus = tempsum_plus*(1 - etaest)*(1 - gamaest[i]);
                else
                    tempsum_plus = tempsum_plus*(1 - gamaest[i]);
            }
            DSMatrix[k][c] = tempsum_plus;
        }
    }
    
    return 0;
}

int KickoutFolds(vector<int> openidx,vector<vector<bool>>Folds,vector<vector<bool>>&FoldsReduce,vector<int>&FoldsReduceIdx)
{
    int FoldsNum = (int)Folds.size();
    for(int i = 0; i<FoldsNum; ++i)
    {
        int count = 0;
        for(vector<int>::iterator it = openidx.begin(); it!=openidx.end(); ++it)
        {
            if(Folds[i][(*it)+1]==0)//+1: folds 0~RNALenFold-1 (3' to 5'), beta 0~RNALenFold-2. beta is start form 2nd of 3' to 5'. The last one in 3' is deleted.
                break;
            else
                ++count;
            
        }
        //if(count == openidx.size())
        if(count>=openidx.size()*0.5)
        {
            FoldsReduce.push_back(Folds[i]);
            FoldsReduceIdx.push_back(i);
        }
    }
    return 0;
}


int ConstructFolds(string FileName, int RNALen, int SampleSize, vector<vector<bool>>&Folds, int&FoldsNum)
{
    ifstream myfile;
    myfile.open (FileName.c_str());
    char * dump = new char[6];
    char * Buffer = new char[RNALen];
    vector<bool> FoldsUnit(RNALen,0);//from 3' to 5'
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
            FoldsUnit[RNALen-1-j] = 1;
        else
            FoldsUnit[RNALen-1-j] = 0;
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
                FoldsUnit[RNALen-1-j] = 1;
            else
                FoldsUnit[RNALen-1-j] = 0;
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
            //  add for text--begin
            /*
            int x = 0;
            for(int pp = 17; pp<24; ++pp)
                x += FoldsUnit[pp];
            if(x==7 & FoldsUnit[24]==0&FoldsUnit[16]==0)
               continue;
            x = 0;
            for(int pp = 16; pp<25; ++pp)
                x += FoldsUnit[pp];
            if(x==9 & FoldsUnit[25]==0&FoldsUnit[15]==0)
                continue;
             */
            // end
            
            Folds.push_back(FoldsUnit);
            ++FoldsNum;
        }
    }
    myfile.close();
    delete []dump;
    delete []Buffer;
    return 0;
}



template<typename T>
double mean(vector<T>in)
{
    int len = (int)in.size();
    double sum = accumulate(in.begin(), in.end(), 0.0);
    return sum/len;
}



