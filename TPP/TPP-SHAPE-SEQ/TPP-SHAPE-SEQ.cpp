//
//  TPP-SHAPE-SEQ.cpp
//  ADD-SHAPE-SEQ
//
//  Created by Hua Li on 5/19/17.
//  Copyright © 2017 Hua Li. All rights reserved.
//

#include "TPP-SHAPE-SEQ.hpp"

using namespace std;
string Folder = "/Users/huali/Desktop/Codes/TPP/";
string RNAname = "TPP80_10";


//TPP
int St = 1;
int Ed = 79;
int RNALenFold = 80;//change
int RNALenNNLS = Ed-St +1+1;//with complete reads



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
    
    //TPP REP2
    {18257, 2521,  2872,   945,  1909,   698,   197,    99,   182,   537,   512,   472,   642,  1269,   635,  1775,  1077,   319,   330,  6388,   998,  2727,
        2141,  2203,   957,   282,   396,   210,   756,   800,   948,  1246,  1161,  2587,  2166, 10351,  2140, 11302,   295,   315,   180,   683,   993,
        6413,  2353,   686,  7886,   451,   314,   241,    98,  1645,  2243,   846,   118,   305,   165,  1470,  4997,   417,   973,   162,   263,   949,
        216,    79 ,   239,   138,   303,   285,   712,   184,    28,    35,   578,   172,   171,   209,  1674, 16648
    };

    


    vector<int>SEQMinusReads=
    {
        2410,1481, 2314,  392, 1567,  286,   48,   25,  184,  531,  199,  143,  312,  405,  353,  164,  479,  309,   49,  583,  159,  642,  390,  251,
        499,  117,  136,  136,  611,   418,  305,  472,  623, 1549,  834, 1535,  471, 1879,   45,  156,   68,  341,  429, 2387,  906,  101, 5959,
        128,  136,   55 ,  10,  668,   204,  216,   50,  115,   69, 1393, 3572,  171,  430,   22,   83,  868,  363,  110,   99,   60,  192,  173,
        280,   51,   15 ,  15,  295,   101,  112,   71,  578, 6807
    };
    
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
    KickoutFolds(opensite,Folds, SEQFoldsReduce,SEQFoldsReduceIdx);//Ignore
    
    
    //Construct Design Matrix in SEQ
    int SEQFoldsReduceNum = (int)SEQFoldsReduce.size();
    vector<vector<double>>SEQDSMatrix(RNALenNNLS,vector<double>(SEQFoldsReduceNum,0));
    ConstructDMSEQ(Ed-St+1, SEQFoldsReduce, SEQetaest, SEQgamaest, SEQDSMatrix);
    
    //Response in EQ
    vector<double>SEQResponse(RNALenNNLS,0);
    //from 3' to 5'
    for(int i = St; i<=Ed; ++i)
        SEQResponse[i-St] = SEQPlusReads[RNALenFold-i];//because plus reads already remove the last nt, so it does not need -1
    SEQResponse[Ed] = SEQPlusReads[0];
    
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
    
    return 0;
    
}

int ConstructDMSEQ(int RNALen, vector<vector<bool>>Folds,double etaest, vector<double>gamaest,vector<vector<double>>&DSMatrix)
{
    int FoldsNum = (int)Folds.size();
    //with complete reads
    for (int k = 0; k < RNALen+1; ++k)//change
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



