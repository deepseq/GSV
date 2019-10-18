/*
 
 ReadsAssembly.cpp
 
 This program is to do:
 1. assemble unmapped reads
 2. summary assembly
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
*/

#include "ReadsAssembly.h"
#include "CmdExecution.h"
#include "CurrentTime.h"
#include "StringOperation.h"
#include "BlastFilter.h"

#include <iostream>
#include <fstream>
#include <regex>
#include <vector>
#include <iomanip>
#include <math.h>

using namespace std;

///////////////////
/* program class */
///////////////////

class contigObj
{
public:
    string contigID;
    int contigLength;
    long int contigCount;
    float contigRPKM;
    
public:
    void getInfo(string subID, int subLen, long int subCount, float subRPKM);
    void addRead();
};

void contigObj::getInfo(string subID, int subLen, long int subCount, float subRPKM)
{
    contigID = subID;
    contigLength = subLen;
    contigCount = subCount;
    contigRPKM = subRPKM;
}

void contigObj::addRead()
{
    contigCount++;
}

///////////////////
/* sub functions */
///////////////////

long int countSamReads(long int cSRNum, string cSRIn, string cSROut)
{
    contigObj cSRContig;
    vector<contigObj> cSRContigList;
    
    regex cSRRegStr("@[A-Z][A-Z]\t");
    
    int cSRi;
    int cSRProIndex = 0;
    string cSRLine;
    vector<string> cSRStr;
    vector<string> cSRTmpID;
    vector<string> cSRTmpLen;
    long int cSRTotal = 0;
    
    ifstream cSRFileIn(cSRIn);
    if(!cSRFileIn){printf("%s \n", "Error! Cannot open file to read! [AssemblyReads]"); exit(1);}
    
    ofstream cSRFileOut(cSROut);
    if(!cSRFileOut){printf("%s \n", "Error! Cannot open file to write! [AssemblyReads]"); exit(1);}
    
    while(1)
    {
        getline(cSRFileIn, cSRLine); // careful! EOF check must be after this cmd as otherwise cannot get right EOF position.
        if(cSRFileIn.eof()) break;
        
        cSRLine.erase(cSRLine.find_last_not_of(" \n\r\t")+1);
        
        if(regex_search(cSRLine, cSRRegStr, regex_constants::match_continuous))
        {
            cSRStr = split(cSRLine, "\t");
            cSRTmpID = split(cSRStr[1], ":");
            cSRTmpLen = split(cSRStr[2], ":");
            cSRContig.getInfo(cSRTmpID[1], str2int(cSRTmpLen[1]), 0, 0.00);
            cSRContigList.push_back(cSRContig);
        }
        else
        {
            cSRStr = split(cSRLine, "\t");
            
            if(cSRStr[2].compare("*") != 0)
            {
                for(cSRi = 0; cSRi < int(cSRContigList.size()); cSRi++)
                {
                    if(cSRStr[2].compare(cSRContigList[cSRi].contigID) == 0)
                    {
                        cSRTotal++;
                        cSRContigList[cSRi].addRead();
                        break;
                    }
                }
            }
        }
        
        // Processing Recording
        cSRProIndex++;
        if(cSRProIndex%10000 == 0)
        {
            cout.flush();
            cout << ">> " << currentDateTime().c_str() << " <<\t------------> Processing Reads:" << cSRProIndex << "\r";
        }
    }
    
    cout.flush();
    cout << ">> " << currentDateTime().c_str() << " <<\t------------> Processing Reads:" << cSRProIndex << "\r" << "\n";
    
    // free vector
    vector<string>().swap(cSRStr);
    vector<string>().swap(cSRTmpID);
    vector<string>().swap(cSRTmpLen);
    
    cSRFileIn.close();
    
    cSRTotal = cSRTotal + cSRNum;
    
    for(cSRi = 0; cSRi < int(cSRContigList.size()); cSRi++)
    {
        cSRFileOut << cSRContigList[cSRi].contigID << "\t" << cSRContigList[cSRi].contigLength << "\t" << cSRContigList[cSRi].contigCount << "\t";
        cSRContigList[cSRi].contigRPKM = (float)cSRContigList[cSRi].contigCount*1000*1000000/(cSRContigList[cSRi].contigLength*cSRTotal);
        cSRFileOut << fixed << setprecision(2) << cSRContigList[cSRi].contigRPKM << "\t";
        if(cSRContigList[cSRi].contigRPKM == 0)
        {
            cSRFileOut << "NA" << "\n";
        }
        else
        {
            cSRFileOut << fixed << setprecision(2) << log10(cSRContigList[cSRi].contigRPKM) << "\n";
        }
    }
    
    cSRFileOut.close();
    
    return cSRTotal;
}

////////////////////
/* main functions */
////////////////////

long int readsAss(string rALen, string rAIns, string rAsdIns, string rARef, string rAID, long int rANum, int rAGap, int rADist, int rAThread)
{
    string readsAssLogStr = rAID + ".log.txt";
    
    long int readsAssTotal = 0;
    string readsAssInTmp;
    string readsAssOutTmp;
    
    string readsFragSize_1 = int2str(str2int(rALen) + str2int(rALen) + str2int(rAIns) - str2int(rAsdIns));
    string readsFragSize_2 = int2str(str2int(rALen) + str2int(rALen) + str2int(rAIns) + str2int(rAsdIns));
    
    string readsAssTmp;
    char *readsAssCmd = new char[1024];
    
    cout << ">> " << currentDateTime().c_str() << " <<\t--------> Assembling Reads\n";
    
    readsAssTmp = "clc_assembler -o " + rAID + "_GSVMiningRes/AssemblyRes/" + rAID + ".contig.tmp.fasta -p fb ss " + readsFragSize_1 + " " + readsFragSize_2 + " -q -i " + rAID + "_GSVMiningRes/UnmappedReads/" + rAID + ".unmapped.read_1.fastq " + rAID + "_GSVMiningRes/UnmappedReads/" + rAID + ".unmapped.read_2.fastq 1>> " + readsAssLogStr;
    strcpy(readsAssCmd, readsAssTmp.c_str());
    cmdExec(readsAssCmd, 0);
    
    string readsAssLine;
    string readsAssFileStr;

    readsAssFileStr = rAID + "_GSVMiningRes/AssemblyRes/" + rAID + ".contig.tmp.fasta";
    ifstream readsAssFileIn(readsAssFileStr);
    if(!readsAssFileIn){printf("%s \n", "Error! Cannot open file to read! [ReadsAssembly]"); exit(1);}
    
    readsAssFileStr = rAID + "_GSVMiningRes/AssemblyRes/" + rAID + ".contig.fasta";
    ofstream readsAssFileOut(readsAssFileStr);
    if(!readsAssFileOut){printf("%s \n", "Error! Cannot open file to write! [ReadsAssembly]"); exit(1);}
    
    string readsAssSeq = "";
    regex readsAssRegStr(">contig_");
    while(1)
    {
        getline(readsAssFileIn, readsAssLine); // careful! EOF check must be after this cmd as otherwise cannot get right EOF position.
        if(readsAssFileIn.eof()) break;
        readsAssLine.erase(readsAssLine.find_last_not_of(" \n\r\t")+1);

        if(regex_search(readsAssLine, readsAssRegStr, regex_constants::match_continuous))
        {
            if(readsAssSeq != "")
            {
                readsAssFileOut << readsAssSeq << "\n";
                readsAssSeq = "";
            }
            readsAssFileOut << readsAssLine << "\n";
        }
        else
        {
            readsAssSeq = readsAssSeq + readsAssLine;
        }
    }
    readsAssFileOut << readsAssSeq << "\n";
    
    readsAssFileOut.close();
    readsAssFileIn.close();
    
    readsAssTmp = "rm -rf " + rAID + "_GSVMiningRes/AssemblyRes/" + rAID + ".contig.tmp.fasta";
    strcpy(readsAssCmd, readsAssTmp.c_str());
    cmdExec(readsAssCmd, 0);
    
    cout << ">> " << currentDateTime().c_str() << " <<\t--------> Blasting Contigs\n";
    
    readsAssTmp = "blastn -db " + rARef + " -query " + rAID + "_GSVMiningRes/AssemblyRes/" + rAID + ".contig.fasta -out " + rAID + "_GSVMiningRes/BlastRes/" + rAID + ".contig.blast.txt -outfmt \"7\"";
    strcpy(readsAssCmd, readsAssTmp.c_str());
    cmdExec(readsAssCmd, 0);
    
    // call function in BlastFilter
    
    GSVBlastFilter(rAID, rAGap, rADist);
    
    cout << ">> " << currentDateTime().c_str() << " <<\t--------> Re-mapping Reads\n";
    
    readsAssTmp = "bwa index " + rAID + "_GSVMiningRes/AssemblyRes/" + rAID + ".GRearr.candidates.fasta 2>> " + readsAssLogStr;
    strcpy(readsAssCmd, readsAssTmp.c_str());
    cmdExec(readsAssCmd, 0);
    
    readsAssTmp = "bwa mem -t " + int2str(rAThread) + " " + rAID + "_GSVMiningRes/AssemblyRes/" + rAID + ".GRearr.candidates.fasta " + rAID + "_GSVMiningRes/UnmappedReads/" + rAID + ".unmapped.read_1.fastq " + rAID + "_GSVMiningRes/UnmappedReads/" + rAID + ".unmapped.read_2.fastq > " + rAID + "_GSVMiningRes/ReMappingRes/" + rAID + ".remap.sam 2>> " + readsAssLogStr;
    strcpy(readsAssCmd, readsAssTmp.c_str());
    cmdExec(readsAssCmd, 0);
    
    readsAssTmp = "awk '{if(($2==99)||($2==147)||($2==83)||($2==163)||($1~/^@/)){print $0}}' " + rAID + "_GSVMiningRes/ReMappingRes/" + rAID + ".remap.sam | grep \"XS:i:0\\|^@\" > " + rAID + "_GSVMiningRes/ReMappingRes/" + rAID + ".remap.unique.sam";
    strcpy(readsAssCmd, readsAssTmp.c_str());
    cmdExec(readsAssCmd, 0);
    
    readsAssTmp = "samtools view -bS " + rAID + "_GSVMiningRes/ReMappingRes/" + rAID + ".remap.unique.sam > " + rAID + "_GSVMiningRes/ReMappingRes/" + rAID + ".remap.unique.bam";
    strcpy(readsAssCmd, readsAssTmp.c_str());
    cmdExec(readsAssCmd, 0);
    
    readsAssTmp = "samtools sort -@ " + int2str(rAThread) + " " + rAID + "_GSVMiningRes/ReMappingRes/" + rAID + ".remap.unique.bam -o " + rAID + "_GSVMiningRes/ReMappingRes/" + rAID + ".remap.unique.sorted.bam";
    strcpy(readsAssCmd, readsAssTmp.c_str());
    cmdExec(readsAssCmd, 0);
    
    readsAssTmp = "samtools index " + rAID + "_GSVMiningRes/ReMappingRes/" + rAID + ".remap.unique.sorted.bam";
    strcpy(readsAssCmd, readsAssTmp.c_str());
    cmdExec(readsAssCmd, 0);
    
    readsAssInTmp = rAID + "_GSVMiningRes/ReMappingRes/" + rAID + ".remap.unique.sam";
    readsAssOutTmp = rAID + "_GSVMiningRes/ReMappingRes/" + rAID + ".contigs.rpkm.txt";
    
    cout << ">> " << currentDateTime().c_str() << " <<\t--------> Re-counting Reads\n";
    
    readsAssTotal = countSamReads(rANum, readsAssInTmp, readsAssOutTmp);
    
    readsAssTmp = "rm -rf " + rAID + "_GSVMiningRes/ReMappingRes/" + rAID + ".remap.sam";
    strcpy(readsAssCmd, readsAssTmp.c_str());
    cmdExec(readsAssCmd, 0);
    
    readsAssTmp = "rm -rf " + rAID + "_GSVMiningRes/ReMappingRes/" + rAID + ".remap.unique.sam";
    strcpy(readsAssCmd, readsAssTmp.c_str());
    cmdExec(readsAssCmd, 0);
    
    readsAssTmp = "rm -rf " + rAID + "_GSVMiningRes/ReMappingRes/" + rAID + ".remap.unique.bam";
    strcpy(readsAssCmd, readsAssTmp.c_str());
    cmdExec(readsAssCmd, 0);
    
    delete[] readsAssCmd;
    
    return readsAssTotal;
}
