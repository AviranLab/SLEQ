//  Created by Hua Li
//  Copyright © 2016 Hua Li. All rights reserved.

#include "DMSTruncation.hpp"

using namespace std;

using namespace std;
///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////         Parameter initialization           ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
string Allele_A = "TGCTGCCATCTCTTTTCTTCTCTATGCGAGGATTTGGACTGGCAGTG";
string Allele_C = "ATCTCTTTTCTTCTCTCTGCGAGGATTTGGACTGGCAGTGAGAATAAGAGACAA";
string CommonSeq = "TTCTTCTCTMTGCGAGGATTTGGACTGGCAGTG";
vector<int>Aidx_AlleleA;//{9,14,17,23,29}
vector<int>Cidx_AlleleA;//{2,5,7,12,24,28}
vector<int>Aidx_AlleleC;//{14,17,23,29}
vector<int>Cidx_AlleleC;//{2,5,7,9,12,24,28}
int RNALenFoldA = 47;
int RNALenFoldC = 54;
int RNALen = 33;
int stOverlapA = 14;
int stOverlapC = 7;
int edOverlapA = 46;
int edOverlapC = 39;
double RhoA = 0.678;
double RhoC = 0.322;


string RNAname = "Allele";
string Folder = "/Users/huali/Desktop/Codes/Human_riboSNitch_alleles/DMSTrunc/"; //!!!!!!!!!!NOTE: Set your own path please.!!!!!!!!!!!!


int main(int argc, const char * argv[]) {
    ofstream myfilew;
    ifstream myfile;
    
    for(int i = stOverlapA; i<=edOverlapA; ++i)
    {
        if(Allele_A[i]=='A')
            Aidx_AlleleA.push_back(i-stOverlapA);
        else if(Allele_A[i]=='C')
            Cidx_AlleleA.push_back(i-stOverlapA);
    }
    for(int i = stOverlapC; i<=edOverlapC; ++i)
    {
        if(Allele_C[i]=='A')
            Aidx_AlleleC.push_back(i-stOverlapC);
        else if(Allele_C[i]=='C')
            Cidx_AlleleC.push_back(i-stOverlapC);
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////           INPUT 1 -- READS            ///////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    int ReadsNum = 0;
    int PatternsNum = 614;//12289;//614;
    vector<Pattern>ReadsAC;
    myfile.open(Folder+"AlleleReadsAC");
    if(myfile)
    {
        for(int i = 0; i<PatternsNum;++i)
        {
            char buffer;
            Pattern tmp;
            for(int j = 0; j<RNALen; ++j)
            {
                myfile>>buffer;
                tmp.pat.push_back(bool(buffer-'0'));
            }
            myfile>>tmp.count;
            ReadsAC.push_back(tmp);
        }
    }
    myfile.close();
    //Reduce mutation patterns into truncation patterns
    vector<int>ReadsSeq(RNALen+1,0);//complete,5'to3'
    for(int i = 0; i<ReadsAC.size(); ++i)
    {
        int j = (int)ReadsAC[i].pat.size()-1;
        for(; j>=0; --j)
            if(ReadsAC[i].pat[j]==1)
            {
                ReadsSeq[j+1] += ReadsAC[i].count;
                break;
            }
        if(j == -1)
            ReadsSeq[0] += ReadsAC[i].count;
    }
    for(int i = 0; i<ReadsSeq.size(); ++i)
        ReadsNum += ReadsSeq[i];
    
    
    
    ////////////////////////////// Downsampling to compare the performance between SEQ and MAP -- Start
    // First calculate the frequencies of all patterns from the data
    /*
    vector<double>Pattern_Freq;
    vector<int>ReadCounts;
    for(int i = 0; i<ReadsSeq.size();++i)
    {
        Pattern_Freq.push_back(double(ReadsSeq[i])/double(ReadsNum));
    }
    //Draw reads according to above frequencies
    int ReadsNum_test = 100;
    vector<double>RandUnifData(ReadsNum_test,0);
    runif(ReadsNum_test,0.0,1.0,RandUnifData);
    sort(RandUnifData.begin(),RandUnifData.end());
    DrawReadsSEQ(RNALen,ReadsNum_test, RandUnifData,Pattern_Freq, ReadCounts);
    //Reput them into ReadsAC
    for(int i = 0; i<ReadsSeq.size(); ++i)
    {
        ReadsSeq[i]= ReadCounts[i];
    }
    ReadsNum = ReadsNum_test;
     */
    ////////////////////////////// Downsampling to compare the performance between SEQ and MAP -- end
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       INPUT 2 -- CANDIDATE STRUCTURES       ///////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //Read dot files of candiate structure of two alleles and remove the duplicates
    int SampleSize = 1000;
    int FoldsNumA = 0;
    int FoldsNumC = 0;
    vector<vector<bool>>FoldsA;//from 3' to 5'
    vector<vector<bool>>FoldsC;
    //Read files & get unique folds & get indices of dominant folds
    string FileNameA = Folder + RNAname+"_A10.txt";
    ConstructFolds(FileNameA, RNALenFoldA, SampleSize, stOverlapA, edOverlapA, FoldsA, FoldsNumA);
    string FileNameC = Folder + RNAname+"_C10.txt";
    ConstructFolds(FileNameC, RNALenFoldC, SampleSize, stOverlapC, edOverlapC, FoldsC, FoldsNumC);
    
    
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////      ETA & BETA ESTIMATION       ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //Order: beta and gama here are from 5' to 3'. Here we do not need reorder beta
    vector<double>SEQbetaest(RNALen,0);//remove the last nucleotide
    for(int i = 0; i<RNALen; ++i)
    {
        //Must have complete reads
        double tempX = (double)accumulate(ReadsSeq.begin(), ReadsSeq.begin()+i+2, 0);
        SEQbetaest[RNALen-i-1] = double(ReadsSeq[i+1])/tempX;
        if(SEQbetaest[RNALen-i-1]<0) SEQbetaest[RNALen-i-1]=0;
    }
    // Manually caluclate EtaA and EtaC
    double EtaA = (SEQbetaest[15]+SEQbetaest[9])/2/RhoC;
    double EtaC = (SEQbetaest[30]/RhoC+(SEQbetaest[27]+SEQbetaest[25]+SEQbetaest[20]+SEQbetaest[8])/RhoA)/5;
    
    
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////      DESIGN MATRIX CONSTRUCTION      ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    int FoldsNum = FoldsNumA+FoldsNumC;
    vector<vector<double>>SEQDSMatrix(RNALen+1,vector<double>(FoldsNum,0));
    ConstructDMSEQ(RNALen, FoldsA, FoldsC, EtaA, EtaC, SEQDSMatrix);
    
   
    //Calculate Y
    vector<double>SEQResponse(RNALen+1,0);//3' to 5', complete
    for(int i = 0; i<RNALen+1; ++i)
        SEQResponse[i] = (double)ReadsSeq[RNALen-i]/(double)ReadsNum;
    
    
    //Write Y
    myfilew.open (Folder+RNAname+"SEQResponse.txt");
    for(int i = 0; i<RNALen+1;++i)
        myfilew<<SEQResponse[i]<<"   ";
    myfilew.close();
    
    
    //Write Design Matrix
    myfilew.open(Folder+RNAname+"SEQDSMatrix.txt");
    for(int i = 0; i<RNALen+1; ++i)
        for(int j = 0; j<FoldsNum; ++j)
            myfilew<<SEQDSMatrix[i][j]<<"    ";
    myfilew.close();
    
    return 0;
    
}
int ConstructDMSEQ(int RNALen, vector<vector<bool>>FoldsA,vector<vector<bool>>FoldsC,double EtaA, double EtaC, vector<vector<double>>&DSMatrix)
{
    int FoldsNumA = (int)FoldsA.size();
    int FoldsNumC = (int)FoldsC.size();
    //with complete reads
    for (int k = 0; k < RNALen+1; ++k)//k starts from 3'
    {
        //FoldsA
        //site k is A/C/UG
        if(k==RNALen)
        {
            for (int c = 0; c < FoldsNumA; ++c)
            {
                float tempsum_plus = 1;
                for (int i = 0; i < k; ++i) //Prob(no adducts on sites from 0 to k-1)
                {
                    if (FoldsA[c][i+(RNALenFoldA-1-edOverlapA)] == 1)
                    {
                        //site i is A/C/UG
                        int flagACbefore = 0;
                        for(vector<int>::iterator itA = Aidx_AlleleA.begin(); itA!=Aidx_AlleleA.end(); ++itA)
                        {
                            if(RNALen-1-*itA == i)
                            {
                                flagACbefore = 1;
                                break;
                            }
                        }
                        for(vector<int>::iterator itC = Cidx_AlleleA.begin(); itC!=Cidx_AlleleA.end() & (flagACbefore == 0); ++itC)
                        {
                            if(RNALen-1-*itC == i)
                            {
                                flagACbefore = 2;
                                break;
                            }
                        }
                        //if UG, tempsum = tempsum * 1;
                        if(flagACbefore == 1)
                            tempsum_plus = tempsum_plus*(1 - EtaA);
                        else if (flagACbefore == 2)
                            tempsum_plus = tempsum_plus*(1 - EtaC);
                    }
                }
                DSMatrix[k][c] = tempsum_plus;
            }
            for(int c = 0; c<FoldsNumC; ++c)
            {
                float tempsum_plus = 1;
                for (int i = 0; i < k; ++i) //Prob(no adducts on sites from 0 to k-1)
                {
                    if (FoldsC[c][i+(RNALenFoldC-1-edOverlapC)] == 1)
                    {
                        int flagACbefore = 0;
                        for(vector<int>::iterator itA = Aidx_AlleleC.begin(); itA!=Aidx_AlleleC.end(); ++itA)
                        {
                            if(RNALen-1-*itA == i)
                            {
                                flagACbefore = 1;
                                break;
                            }
                        }
                        for(vector<int>::iterator itC = Cidx_AlleleC.begin(); itC!=Cidx_AlleleC.end()&(flagACbefore==0); ++itC)
                        {
                            if(RNALen-1-*itC == i)
                            {
                                flagACbefore = 2;
                                break;
                            }
                        }
                        if(flagACbefore == 1)
                            tempsum_plus = tempsum_plus*(1 - EtaA);
                        else if (flagACbefore == 2)
                            tempsum_plus = tempsum_plus*(1 - EtaC);
                    }
                }
                DSMatrix[k][c+FoldsNumA] = tempsum_plus;
            }
            continue;
        }
        int flagAC = 0;
        for(vector<int>::iterator itA = Aidx_AlleleA.begin(); itA!=Aidx_AlleleA.end(); ++itA)
        {
            if(RNALen-1-*itA == k)
            {
                flagAC = 1;
                break;
            }
        }
        for(vector<int>::iterator itC = Cidx_AlleleA.begin(); itC!=Cidx_AlleleA.end() & (flagAC == 0); ++itC)
        {
            if(RNALen-1-*itC == k)
            {
                flagAC = 2;
                break;
            }
        }
        if(flagAC == 0)//UG,if in UG, mutation is impossible since we processed it.
        {//Prob of forming an adduct on UG sites is 0;//already initialized 0
            //for(int c = 0; c<FoldsNumA+FoldsNumC; ++c)
              //  DSMatrix[k][c] = 0;
            continue;
        }
        else
        {
            for (int c = 0; c < FoldsNumA; ++c)
            {
                float tempsum_plus = 1;
                if(k < RNALen) //Prob(adduct on site k)
                {
                    if(FoldsA[c][k+(RNALenFoldA-1-edOverlapA)] == 1)
                    {
                        if(flagAC == 1)
                            tempsum_plus = EtaA;
                        else if (flagAC == 2)
                            tempsum_plus = EtaC;
                    }
                    else
                    {
                        //tempsum_plus = 0;//initialization is 0
                        continue;
                    }
                }
                for (int i = 0; i < k; ++i) //Prob(no adducts on sites from 0 to k-1)
                {
                    if (FoldsA[c][i+(RNALenFoldA-1-edOverlapA)] == 1)
                    {
                        //site i is A/C/UG
                        int flagACbefore = 0;
                        for(vector<int>::iterator itA = Aidx_AlleleA.begin(); itA!=Aidx_AlleleA.end(); ++itA)
                        {
                            if(RNALen-1-*itA == i)
                            {
                                flagACbefore = 1;
                                break;
                            }
                        }
                        for(vector<int>::iterator itC = Cidx_AlleleA.begin(); itC!=Cidx_AlleleA.end() & (flagACbefore == 0); ++itC)
                        {
                            if(RNALen-1-*itC == i)
                            {
                                flagACbefore = 2;
                                break;
                            }
                        }
                        //if UG, tempsum = tempsum * 1;
                        if(flagACbefore == 1)
                            tempsum_plus = tempsum_plus*(1 - EtaA);
                        else if (flagACbefore == 2)
                            tempsum_plus = tempsum_plus*(1 - EtaC);
                    }
                }
                DSMatrix[k][c] = tempsum_plus;
            }
        }

        //FoldsC
        flagAC = 0;
        for(vector<int>::iterator itA = Aidx_AlleleC.begin(); itA!=Aidx_AlleleC.end(); ++itA)
        {
            if(RNALen-1-*itA == k)
            {
                flagAC = 1;
                break;
            }
        }
        for(vector<int>::iterator itC = Cidx_AlleleC.begin(); itC!=Cidx_AlleleC.end()& (flagAC==0); ++itC)
        {
            if(RNALen-1-*itC == k)
            {
                flagAC = 2;
                break;
            }
        }
        if(flagAC != 0)
        {
            for(int c = 0; c<FoldsNumC; ++c)
            {
                float tempsum_plus = 1;
                if(k < RNALen) //Prob(adduct on site k)
                {
                    if(FoldsC[c][k+(RNALenFoldC-1-edOverlapC)] == 1)
                    {
                        if(flagAC == 1)
                            tempsum_plus = EtaA;
                        else if (flagAC == 2)
                            tempsum_plus = EtaC;
                    }
                    else
                        continue;
                        //tempsum_plus = 0;
                }
                for (int i = 0; i < k; ++i) //Prob(no adducts on sites from 0 to k-1)
                {
                    if (FoldsC[c][i+(RNALenFoldC-1-edOverlapC)] == 1)
                    {
                        int flagACbefore = 0;
                        for(vector<int>::iterator itA = Aidx_AlleleC.begin(); itA!=Aidx_AlleleC.end(); ++itA)
                        {
                            if(RNALen-1-*itA == i)
                            {
                                flagACbefore = 1;
                                break;
                            }
                        }
                        for(vector<int>::iterator itC = Cidx_AlleleC.begin(); itC!=Cidx_AlleleC.end()&(flagACbefore==0); ++itC)
                        {
                            if(RNALen-1-*itC == i)
                            {
                                flagACbefore = 2;
                                break;
                            }
                        }
                        if(flagACbefore == 1)
                            tempsum_plus = tempsum_plus*(1 - EtaA);
                        else if (flagACbefore == 2)
                            tempsum_plus = tempsum_plus*(1 - EtaC);
                    }
                }
                DSMatrix[k][c+FoldsNumA] = tempsum_plus;
            }
        }
    }
    
    return 0;
}


int ConstructFolds(string FileName, int RNALen, int SampleSize, int st, int ed, vector<vector<bool>>&Folds, int&FoldsNum)
{
    ifstream myfile;
    myfile.open (FileName.c_str());
    char * dump = new char[6];
    char * Buffer = new char[RNALen];
    vector<bool> FoldsUnit(RNALen,0);//from 3' to 5'
    FoldsNum = 0;
    
    int dumpnum = 0;
    //Read Files, convert to 0,1 form, remain unique folds
    for(int i = 0; i<20; ++i)
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
            /*if(Folds[k] == FoldsUnit)
             {
             flag = 1;
             break;
             }*/
            int q = RNALen-ed-1;
            for(; q<=RNALen-st-1; ++q)
            {
                if(Folds[k][q]!=FoldsUnit[q])
                    break;
            }
            if(q == RNALen-st)
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



int DrawReadsSEQ(int RNALen,int ReadsNum, vector<double>RandUnifData,vector<double>&ProbPlus, vector<int>& PlusReads)
{
    
    vector<double>ProbPlusCumsum(RNALen+1,0);
    partial_sum(ProbPlus.begin(), ProbPlus.end(), ProbPlusCumsum.begin());
    
    int PlusBegin = 0;
    for(int i = 0; i<RNALen+1; ++i)
    {
        for(int j = PlusBegin; j<ReadsNum; ++j)
        {
            if(RandUnifData[j]>ProbPlusCumsum[i])
            {
                PlusReads.push_back(j - PlusBegin);
                PlusBegin = j;
                break;
            }
            if(j == ReadsNum-1)
            {
                PlusReads.push_back(ReadsNum - PlusBegin);
                PlusBegin = ReadsNum;
            }
        }
    }
    
    for(int i = (int)PlusReads.size(); i<RNALen+1; ++i)
        PlusReads.push_back(0);

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