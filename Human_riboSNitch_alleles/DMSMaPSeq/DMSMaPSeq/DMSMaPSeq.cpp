//  Created by Hua Li
//  Copyright © 2016 Hua Li. All rights reserved.


#include "DMSMaPSeq.hpp"

using namespace std;
///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////         Parameter initialization           ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
string Allele_A = "TGCTGCCATCTCTTTTCTTCTCTATGCGAGGATTTGGACTGGCAGTG"; //sequence of allele A
string Allele_C = "ATCTCTTTTCTTCTCTCTGCGAGGATTTGGACTGGCAGTGAGAATAAGAGACAA"; //sequence of allele C
string CommonSeq = "TTCTTCTCTMTGCGAGGATTTGGACTGGCAGTG"; // common sequence of two alleles
vector<int>Aidx_AlleleA;//{9,14,17,23,29} // Positions of adenine in allele A
vector<int>Cidx_AlleleA;//{2,5,7,12,24,28} // Positions of cytosine in allele A
vector<int>Aidx_AlleleC;//{14,17,23,29} // Positions of adenine in allele C
vector<int>Cidx_AlleleC;//{2,5,7,9,12,24,28} // Positions of cytosine in allele C
int RNALenFoldA = 47; //RNA length of allele A
int RNALenFoldC = 54; //RNA length of allele C
int RNALen = 33; // common length of two alleles
int stOverlapA = 14; // start position of common sequence in allele A
int stOverlapC = 7; // start position of common sequence in allele C
int edOverlapA = 46; // end position of common sequence in allele A
int edOverlapC = 39; // end position of common sequence in allele C
double RhoA = 0.678; //relative abundance of allele A
double RhoC = 0.322; //relative abundance of allele C

string RNAname = "Allele"; // RNA name
string Folder = "/Users/huali/Desktop/Codes/Human_riboSNitch_alleles/DMSMaPSeq/"; //!!!!!!!!!!NOTE: Set your own path please.!!!!!!!!!!!!


int main(int argc, const char * argv[]) {
    
    ofstream myfilew;
    ifstream myfile;
    
    // Select positions of adenine and cytosine in alleles A and C
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
    
    //find indices of UG sites
    vector<int>UGidx;
    for(int i = 0; i<RNALen; ++i)
    {
        if(CommonSeq[i]=='G'||CommonSeq[i]=='T')
            UGidx.push_back(i);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////           INPUT 1 -- READS            ///////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //Bin reads into patterns
    //myfile.open (Folder + "SRR3929636.sam"); // open DMS-MaPSeq file
    //if(!myfile)
    //{
    //    cout<<"Failed to read the file!" <<endl;
    //    return 1;
    //}
    vector<Pattern>Reads;
    string s;
    int ReadsNum = 0;
    int CountA_Tol = 0; // Number of reads aligned to allele A
    int CountC_Tol = 0; // Number of reads aligned to allele C
    int CountA_F = 0; // Number of complete reads aligned to allele A
    int CountC_F = 0; // Number of complete reads aligned to allele C
    int ReadsNumTol = 102466345;
    int stSAM = 12;
    /*
    for(int ridx = 0; ridx<ReadsNumTol; ++ridx)//read iteration
    {
        sam CurrentRead;
        CurrentRead.pos = 0; //initialize
        vector<int> cigarnum; //save numbers in cigar
        vector<char> cigarchar; // save characters in cigar
        vector<int>cigar_cutidx; // save indices of cutoff
        vector<bool>CurrentPat(RNALen,0); //save pattern of current read
        
        int flag = 0;
        int flag2 = 0;
        int flag3 = 0;
        int flag4 = 0;
        int flag5 = 0;
        int flag6 = 0;
        int flag9 = 0;
        int flag10 = 0;
        int breakflagAlgn = 0;
        int breakflagSM = 0;
        int breakflagFull = 0;
        
        int stPos = 0; //start index of postion in the string
        int stCigar = 0; //start index of cigar in the string
        int stSeq = 0; //start index of seq in the string
        int stMatch = 0;//start index of matchx
        
        //Extract information of a read from sam file
        getline(myfile,s); //get one line of sam file, it contains all information of one read
        for(int sidx = stSAM; sidx<s.length(); ++sidx) //sidx: string index
        {
            if(s[sidx]=='\t') // find cutoff points
            {
                ++flag;
            }
            if(flag == 2 & flag2 == 0)// aligned information
            {
                if(s[sidx-1]-'0'!=0) //if not aligend, stop and go to next read;
                {
                    breakflagAlgn = 1;
                    break;
                }
                flag2 = 1;// avoid to change the result
            }
            if(flag == 3 & flag3 ==0)
            {
                stPos = sidx+1;// start index of postion in the string
                if(s[sidx-1]-'C'==0) //aligned to Allele_C
                {
                    ++CountC_Tol;
                    CurrentRead.ref = 'C';
                }
                else if (s[sidx-1]-'A'==0) //aligend to Allele_A
                {
                    ++CountA_Tol;
                    CurrentRead.ref = 'A';
                }
                else
                {cout<<"No Reference!"<<endl;
                    return 1;}
                flag3 = 1;
            }
            if(flag == 4 && flag4 == 0)
            {
                int l = sidx-stPos;//length of postion number
                for(int k = stPos; k<sidx; ++k)
                {
                    CurrentRead.pos += pow(10,l-1)*(int)(s[k]-'0');
                    --l;
                }
                flag4 = 1;
            }
            if(flag == 5 & flag5 == 0)
            {
                stCigar = sidx+1;//start index of cigar in the string
                flag5 = 1;
            }
            if(flag == 6 & flag6 == 0)
            {
                for(int k = stCigar; k<sidx; ++k)
                    CurrentRead.cigar.push_back(s[k]);//cigar
                flag6 = 1;
            }
            
            if(flag == 9 & flag9 == 0)
            {
                stSeq = sidx+1;//start index of seq in the string
                flag9 = 1;
            }
            if(flag == 10 & flag10 == 0)
            {
                for(int k = stSeq; k<sidx; ++k)
                    CurrentRead.seq.push_back(s[k]);//seq
                flag10 = 1;
                break;//stop string iteration
            }
        }
        if(breakflagAlgn == 1)
            continue; //if not aligned, continue;
        if(CurrentRead.ref-'A'!=0&CurrentRead.ref-'C'!=0)
        {
            continue;
        }
        
        //Cut cigar
        cigar_cutidx.push_back(-2);
        for(int cidx = 0; cidx<CurrentRead.cigar.size(); ++cidx)//cidx: cigar index
        {
            if(CurrentRead.cigar[cidx]<'0' || CurrentRead.cigar[cidx]>'9')//Not a number
            {
                if(CurrentRead.cigar[cidx]-'M'!=0 & CurrentRead.cigar[cidx]-'S'!=0)//contains character besides 'M' and 'S'
                {
                    breakflagSM = 1;//break current cigar cut and continue read iteration
                    break;
                }
                cigar_cutidx.push_back(cidx-1);//if it is'M' or 'S', save the index of end of the number
                cigarchar.push_back(CurrentRead.cigar[cidx]);// and save the character
            }
        }
        if(breakflagSM == 1)// contains character besides 'M' and 'S', continue read iteration
            continue;
        int flagM = 0;
        if(CurrentRead.ref-'A'==0)//aligned to A
        {
            if(CurrentRead.pos-1>stOverlapA)//start postion is larger than start postin of common seq
                continue;
            else
            {
                //construct cigar num
                for(int i = 1; i<cigar_cutidx.size(); ++i)//i: iterate all cut indices
                {
                    int l = cigar_cutidx[i]-cigar_cutidx[i-1]-1;//length of number
                    int tempnum = 0;
                    for(int k = cigar_cutidx[i-1]+2; k<=cigar_cutidx[i]; ++k)
                    {
                        tempnum += pow(10,l-1)*(int)(CurrentRead.cigar[k]-'0');
                        --l;
                    }
                    cigarnum.push_back(tempnum);
                }
                
                //Compare seq
                if(cigarnum.size()!=cigarchar.size())
                {
                    cout<<"Mismatch of cigar chars and nums!"<<endl;
                    return 1;
                }
                else
                {
                    for(int i = 0; i<cigarchar.size()&flagM==0; ++i)
                    {
                        if(cigarchar[i]-'M' == 0)
                        {
                            if(cigarnum[i]+CurrentRead.pos-1<edOverlapA+1)//st(real) + readlength -1 = ed(real) = ed(c)+1, end position is less than that of overlap
                            {
                                breakflagFull = 1;//break current iteration
                                break;
                            }
                            flagM = 1;//not allow multiple 'M'
                        }
                        else stMatch += cigarnum[i];
                    }
                    if(breakflagFull==1)
                        continue;//continue read iteration
                    ++CountA_F;
                    for(int i = stOverlapA; i<=edOverlapA; ++i)//compare read seq and reference
                    {
                        if(CurrentRead.seq[i+2-CurrentRead.pos+stMatch-1]!=Allele_A[i])//l = ed - st +1; the other +1 is for i is c end; st+l - 1 = ed
                        {
                            CurrentPat[i-stOverlapA] = 1;
                        }
                    }
                }
            }
        }
        if(CurrentRead.ref-'C'==0)
        {
            if(CurrentRead.pos-1>stOverlapC)
                continue;
            else
            {
                //cigar num
                for(int i = 1; i<cigar_cutidx.size(); ++i)
                {
                    int l = cigar_cutidx[i]-cigar_cutidx[i-1]-1;
                    int tempnum = 0;
                    for(int k = cigar_cutidx[i-1]+2; k<=cigar_cutidx[i]; ++k)
                    {
                        tempnum += pow(10,l-1)*(int)(CurrentRead.cigar[k]-'0');
                        --l;
                    }
                    cigarnum.push_back(tempnum);
                }
                
                //Compare seq
                if(cigarnum.size()!=cigarchar.size())
                {
                    cout<<"Mismatch of cigar chars and nums!"<<endl;
                    return 1;
                }
                else
                {
                    for(int i = 0; i<cigarchar.size()&flagM==0; ++i)
                    {
                        if(cigarchar[i]-'M' == 0)
                        {
                            if(cigarnum[i]+CurrentRead.pos-1<edOverlapC+1)
                            {
                                breakflagFull = 1;
                                break;
                            }
                            flagM=1;
                        }
                        else stMatch += cigarnum[i];
                    }
                    if(breakflagFull==1)
                        continue;
                    ++CountC_F;
                    for(int i = stOverlapC; i<=edOverlapC; ++i)
                    {
                        if(CurrentRead.seq[i+2-CurrentRead.pos+stMatch-1]!=Allele_C[i])//i+2-CurrentRead.pos+stMatch-1
                        {
                            CurrentPat[i-stOverlapC] = 1;
                        }
                    }
                }
            }
        }
        
        int iexist = 0;
        for(; iexist<Reads.size(); ++iexist)//check if current read pattern exist
        {
            if(Reads[iexist].pat==CurrentPat)//if exist, only update the number of reads of this pattern
            {
                Reads[iexist].count = Reads[iexist].count+1;//check initilize
                break;
            }
        }
        if(iexist == Reads.size())//new pattern.
        {
            Pattern newpat;
            newpat.pat = CurrentPat;
            newpat.count = 1;
            Reads.push_back(newpat);
        }
    }
    //myfile.close();
    //Ignore mutations at U and G sites and re-bin reads into patterns
    vector<Pattern>ReadsAC;
    
    for(int i = 0; i<Reads.size(); ++i)
    {
        ReadsNum += Reads[i].count;
        for(vector<int>::iterator it = UGidx.begin(); it!=UGidx.end(); ++it)
        {
            if(Reads[i].pat[*it]==1)
                Reads[i].pat[*it] = 0;
        }
        int iexist = 0;
        for(; iexist<ReadsAC.size(); ++iexist)//check if current read pattern exist
        {
            if(ReadsAC[iexist].pat==Reads[i].pat)//if exist, only update the number of reads of this pattern
            {
                ReadsAC[iexist].count = ReadsAC[iexist].count+Reads[i].count;//check initilize
                break;
            }
        }
        if(iexist == ReadsAC.size())//new pattern. //check - 1 or
        {
            for(int k = 0; k<RNALen; ++k)
                if(Reads[i].pat[k]==1)
                    Reads[i].muidx.push_back(k);
            ReadsAC.push_back(Reads[i]);
        }
    }
    int PatternsNum = (int)ReadsAC.size();
    vector<double>Count(10,0);//mu = 0,1,2,3,>=4
    for(int i = 0; i<PatternsNum; ++i)
    {
        int mu = (int)ReadsAC[i].muidx.size();
        Count[mu] += ReadsAC[i].count;
    }
    
    myfilew.open(Folder+RNAname+"ReadsAC");
    for(int i = 0; i<PatternsNum; ++i)
    {
        for(int j = 0; j<ReadsAC[i].pat.size(); ++j)
            myfilew<<ReadsAC[i].pat[j];
        myfilew<<"  "<<ReadsAC[i].count<<endl;
    }
    myfilew.close();
*/
    
    //test_begin
    
    vector<Pattern>ReadsAC;
    myfile.open (Folder+RNAname+"ReadsAC");
    int PatternsNum = 614;
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
            ReadsNum += tmp.count;
            ReadsAC.push_back(tmp);
        }
    }
    myfile.close();
    ////////////////////////////// Downsampling to compare the performance between SEQ and MAP -- Start
    // First calculate the frequencies of all patterns from the data
    /*
    vector<double>Pattern_Freq;
    vector<int>ReadCounts;
    for(int i = 0; i<ReadsAC.size();++i)
    {
        Pattern_Freq.push_back(double(ReadsAC[i].count)/double(ReadsNum));
    }
    //Draw reads according to above frequencies
    int ReadsNum_test = 1000;
    vector<double>RandUnifData(ReadsNum_test,0);
    runif(ReadsNum_test,0.0,1.0,RandUnifData);
    sort(RandUnifData.begin(),RandUnifData.end());
    DrawReadsMAP((int)ReadsAC.size(),ReadsNum_test, RandUnifData,Pattern_Freq, ReadCounts);
    //Reput them into ReadsAC
    for(int i = 0; i<ReadsAC.size(); ++i)
    {
        ReadsAC[i].count = ReadCounts[i];
    }
    ReadsNum = ReadsNum_test;
    ////////////////////////////// Downsampling to compare the performance between SEQ and MAP -- end
     */
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       INPUT 2 -- CANDIDATE STRUCTURES       ///////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //Read dot files of candiate structure of two alleles and remove the duplicates
    int SampleSize = 1000;
    int FoldsNumA = 0;
    int FoldsNumC = 0;
    vector<vector<bool>>FoldsA; //structure candidate set of allele A. from 3' to 5'
    vector<vector<bool>>FoldsC; //structure candidate set of allele C.
    string FileNameA = Folder + RNAname+"_A10.txt";
    ConstructFolds(FileNameA, RNALenFoldA, SampleSize, stOverlapA, edOverlapA, FoldsA, FoldsNumA);
    string FileNameC = Folder + RNAname+"_C10.txt";
    ConstructFolds(FileNameC, RNALenFoldC, SampleSize, stOverlapC, edOverlapC, FoldsC, FoldsNumC);
    
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////      ETA & BETA ESTIMATION       ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //Order: beta and gama here are from 5' to 3'. Here we do not need reorder beta
    vector<double>MAPbetaest(RNALen,0);//remove the last nucleotide
    vector<int>pJc(RNALen,0);
    for(int p = 0; p<ReadsAC.size(); ++p)
    {
        for(int s = 0; s<RNALen; ++s)
        {
            if(ReadsAC[p].pat[s]==1)
            {
                pJc[s] += ReadsAC[p].count;
            }
        }
        for(int i = 0; i<RNALen; ++i)
        {
            MAPbetaest[i] = double(pJc[i])/ReadsNum;
        }
    }
    // Manually caluclate EtaA and EtaC
    double EtaA = (MAPbetaest[17]+MAPbetaest[23])/2/RhoC;
    double EtaC = (MAPbetaest[2]/RhoC+(MAPbetaest[5]+MAPbetaest[7]+MAPbetaest[12]+MAPbetaest[24])/RhoA)/5;
    //EtaA = 0.08;
    //EtaC = 0.056;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////      DESIGN MATRIX CONSTRUCTION      ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    int FoldsNum = FoldsNumA+FoldsNumC;
    vector<vector<float>>MAPDSMatrix(PatternsNum,vector<float>(FoldsNum,0));
    ConstructDMMAP(RNALen, FoldsA, FoldsC, ReadsAC, EtaA, EtaC, MAPDSMatrix);
    

    // Calculate Y
    vector<double>MAPResponse(PatternsNum,0);
    for(int i = 0; i<PatternsNum; ++i)
    {
        MAPResponse[i] = (double)ReadsAC[i].count/(double)ReadsNum;
    }

    //Write Y
    myfilew.open (Folder+RNAname+"MAPResponse.txt");
    for(int i = 0; i<PatternsNum;++i)
        myfilew<<MAPResponse[i]<<"   ";
    myfilew.close();

    
    //Write Design Matrix
    myfilew.open(Folder+RNAname+"MAPDSMatrix.txt");
    for(int i = 0; i<PatternsNum; ++i)
        for(int j = 0; j<FoldsNum; ++j)
            myfilew<<MAPDSMatrix[i][j]<<"    ";
    myfilew.close();
   
    myfilew.open (Folder+RNAname+"SEQProfile.txt");
    for(int i = 0; i<RNALen;++i)
        myfilew<<MAPbetaest[i]<<"   ";
    myfilew.close();
    
    return 0;
}


int ConstructDMMAP(int RNALen, vector<vector<bool>>FoldsA, vector<vector<bool>>FoldsC,vector<Pattern>ReadsAC, double EtaA, double EtaC, vector<vector<float>>&DSMatrix)
{
    int FoldsNumA = (int)FoldsA.size();
    int FoldsNumC = (int)FoldsC.size();
    int PatternsNum = (int)ReadsAC.size();
    for (int p = 0; p < PatternsNum; ++p)
    {
        
        //FoldsA
        for(int c = 0; c<FoldsNumA; ++c)
        {
            float tempsum_plus = 1;
            for(vector<int>::iterator itA = Aidx_AlleleA.begin(); itA!=Aidx_AlleleA.end(); ++itA)
            {
                if(ReadsAC[p].pat[*itA]==0)//no mutation
                {
                    if(FoldsA[c][*itA+stOverlapA] == 1)//unpaired
                        tempsum_plus = tempsum_plus*(1-EtaA);
                }
                else //mutation
                {
                    if(FoldsA[c][*itA+stOverlapA] == 1)//unpaired
                        tempsum_plus = tempsum_plus*EtaA;
                    else
                        tempsum_plus = tempsum_plus*0;
                }
            }
            for(vector<int>::iterator itC = Cidx_AlleleA.begin(); itC!=Cidx_AlleleA.end(); ++itC)
            {
                if(ReadsAC[p].pat[*itC]==0)//no mutation
                {
                    if(FoldsA[c][*itC+stOverlapA] == 1)//unpaired
                        tempsum_plus = tempsum_plus*(1-EtaC);
                }
                else //mutation
                {
                    if(FoldsA[c][*itC+stOverlapA] == 1)//unpaired
                        tempsum_plus = tempsum_plus*EtaC;
                    else
                        tempsum_plus = tempsum_plus*0;
                }
            }
            DSMatrix[p][c] = tempsum_plus;
        }
        
        //FoldsC
        for(int c = 0; c<FoldsNumC; ++c)
        {
            float tempsum_plus = 1;
            for(vector<int>::iterator itA = Aidx_AlleleC.begin(); itA!=Aidx_AlleleC.end(); ++itA)
            {
                if(ReadsAC[p].pat[*itA]==0)//no mutation
                {
                    if(FoldsC[c][*itA+stOverlapC] == 1)//unpaired
                        tempsum_plus = tempsum_plus*(1-EtaA); //unpaired&no mutaiton&A
                }
                else //mutation
                {
                    if(FoldsC[c][*itA+stOverlapC] == 1)//unpaired
                        tempsum_plus = tempsum_plus*EtaA;
                    else
                        tempsum_plus = tempsum_plus*0;
                }
            }
            for(vector<int>::iterator itC = Cidx_AlleleC.begin(); itC!=Cidx_AlleleC.end(); ++itC)
            {
                if(ReadsAC[p].pat[*itC]==0)//no mutation
                {
                    if(FoldsC[c][*itC+stOverlapC] == 1)//unpaired
                        tempsum_plus = tempsum_plus*(1-EtaC);
                }
                else //mutation
                {
                    if(FoldsC[c][*itC+stOverlapC] == 1)//unpaired
                        tempsum_plus = tempsum_plus*EtaC;
                    else
                        tempsum_plus = tempsum_plus*0;
                }
            }
            DSMatrix[p][c+FoldsNumA] = tempsum_plus;
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
    //vector<bool> FoldsUnit(RNALen,0);//from 3' to 5'
    vector<bool> FoldsUnit(RNALen,0);//from 5' to 3'
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
    /*
    for(int j = 0; j<RNALen; ++j)
    {
        if(Buffer[j] == '.')
            FoldsUnit[RNALen-1-j] = 1;
        else
            FoldsUnit[RNALen-1-j] = 0;
    }
    */
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
        /*
        for(int j = 0; j<RNALen; ++j)
        {
            if(Buffer[j] == '.')
                FoldsUnit[RNALen-1-j] = 1;
            else
                FoldsUnit[RNALen-1-j] = 0;
        }
         */
        for(int j = 0; j<RNALen; ++j)
        {
            if(Buffer[j] == '.')
                FoldsUnit[j] = 1;
            else
                FoldsUnit[j] = 0;
        }
        /*
        for(int k = 0; k<FoldsNum; ++k)
        {
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
         */
        for(int k = 0; k<FoldsNum; ++k)
        {
            int q = st;
            for(; q<=ed; ++q)
            {
                if(Folds[k][q]!=FoldsUnit[q])
                    break;
            }
            if(q == ed+1)
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
    
    
    return 0;
}


int DrawReadsMAP(int PatternsNum,int ReadsNum, vector<double>RandUnifData,vector<double>&ProbPlus, vector<int>& PlusReads)
{
    
    vector<double>ProbPlusCumsum(PatternsNum,0);
    partial_sum(ProbPlus.begin(), ProbPlus.end(), ProbPlusCumsum.begin());
    
    int PlusBegin = 0;
    for(int i = 0; i<PatternsNum; ++i)
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
    for (int i = (int)PlusReads.size(); i<PatternsNum; ++i)
    {
        PlusReads.push_back(0);
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