/*
 
 CoverageFilter.cpp
 
 This program is to do:
 1. fitler false positive results
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
*/

#include "CoverageFilter.h"
#include "CurrentTime.h"
#include "CmdExecution.h"
#include "StringOperation.h"
#include "ProcessTool.h"

#include <iostream>
#include <fstream>
#include <regex>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <math.h>

using namespace std;

int getFlankCounts(string gFCBam, string gFCID, long int gFCStart, long int gFCEnd, int gFCWin, long int gFCTotal)
{
    int gFCNum;
    string gFCNumStr;
    
    string gFCCmdTmp;
    char *gFCCmd = new char[1024];
    
    gFCCmdTmp = "samtools view -c " + gFCBam + " " + gFCID + ":" + longint2str(gFCStart) + "-" + longint2str(gFCEnd);
    strcpy(gFCCmd, gFCCmdTmp.c_str());
    gFCNumStr = cmdGetOut(gFCCmd, 0);
    gFCNum = str2int(gFCNumStr);
    // potential problem, float is enough for the small or large number?
    gFCNum = gFCNum/(float(gFCWin)/1000)/(float(gFCTotal)/1000000);
    
    return gFCNum;
}

pair<string, string> flankCalCov(string fCCBwa, string fCCContig, int fCCGroup, string fCCContigID, string fCCChrID, long int fCCContigStart, long int fCCContigEnd, long int fCCChrStart, long int fCCChrEnd, string fCCOri, long int fCCBWATotal, long int fCCContigTotal, int fCCWinSize)
{
    string fCCResA;
    string fCCResB;
    
    long int fCCConStartPos;
    long int fCCConEndPos;
    long int fCCBwaStartPos;
    long int fCCBwaEndPos;
    
    int fCCContigCountIn;
    int fCCContigCountOut;
    int fCCBwaCountIn;
    int fCCBwaCountOut;
    
    string fCCContigTmp;
    string fCCBwaTmp;
    
    if(fCCGroup == 0)
    {
        /* contig coverage */
        fCCConStartPos = fCCContigEnd - fCCWinSize;
        fCCConEndPos = fCCContigEnd + fCCWinSize;
        if(fCCConStartPos < 0) fCCConStartPos = 0;
        
        fCCContigCountIn = getFlankCounts(fCCContig, fCCContigID, fCCConStartPos, fCCContigEnd, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        fCCContigCountOut = getFlankCounts(fCCContig, fCCContigID, fCCContigEnd, fCCConEndPos, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        
        /* genome coverage */
        fCCBwaStartPos = fCCChrEnd - fCCWinSize;
        fCCBwaEndPos = fCCChrEnd + fCCWinSize;
        if(fCCBwaStartPos < 0) fCCBwaStartPos = 0;
        
        if(fCCOri == "forward")
        {
            fCCBwaCountIn = getFlankCounts(fCCBwa, fCCChrID, fCCBwaStartPos, fCCChrEnd, fCCWinSize, (fCCBWATotal+fCCContigTotal));
            fCCBwaCountOut = getFlankCounts(fCCBwa, fCCChrID, fCCChrEnd, fCCBwaEndPos, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        }
        else
        {
            fCCBwaCountIn = getFlankCounts(fCCBwa, fCCChrID, fCCChrEnd, fCCBwaEndPos, fCCWinSize, (fCCBWATotal+fCCContigTotal));
            fCCBwaCountOut = getFlankCounts(fCCBwa, fCCChrID, fCCBwaStartPos, fCCChrEnd, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        }
        
        fCCResA = int2str(fCCContigCountIn) + "--" + int2str(fCCContigCountOut);
        fCCResB = int2str(fCCBwaCountIn) + "--" + int2str(fCCBwaCountOut);
        
    }
    else if (fCCGroup == -1)
    {
        fCCConStartPos = fCCContigStart - fCCWinSize;
        fCCConEndPos = fCCContigStart + fCCWinSize;
        if(fCCConStartPos < 0) fCCConStartPos = 0;
        
        fCCContigCountIn = getFlankCounts(fCCContig, fCCContigID, fCCContigStart, fCCConEndPos, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        fCCContigCountOut = getFlankCounts(fCCContig, fCCContigID, fCCConStartPos, fCCContigStart, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        
        fCCBwaStartPos = fCCChrStart - fCCWinSize;
        fCCBwaEndPos = fCCChrStart + fCCWinSize;
        if(fCCBwaStartPos < 0) fCCBwaStartPos = 0;
        
        if(fCCOri == "forward")
        {
            fCCBwaCountIn = getFlankCounts(fCCBwa, fCCChrID, fCCChrStart, fCCBwaEndPos, fCCWinSize, (fCCBWATotal+fCCContigTotal));
            fCCBwaCountOut = getFlankCounts(fCCBwa, fCCChrID, fCCBwaStartPos, fCCChrStart, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        }
        else
        {
            fCCBwaCountIn = getFlankCounts(fCCBwa, fCCChrID, fCCBwaStartPos, fCCChrStart, fCCWinSize, (fCCBWATotal+fCCContigTotal));
            fCCBwaCountOut = getFlankCounts(fCCBwa, fCCChrID, fCCChrStart, fCCBwaEndPos, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        }
        
        fCCResA = int2str(fCCContigCountIn) + "--" + int2str(fCCContigCountOut);
        fCCResB = int2str(fCCBwaCountIn) + "--" + int2str(fCCBwaCountOut);
    }
    else
    {
        /* contig coverage */
        fCCConStartPos = fCCContigStart - fCCWinSize;
        fCCConEndPos = fCCContigStart + fCCWinSize;
        if(fCCConStartPos < 0) fCCConStartPos = 0;
        
        fCCContigCountIn = getFlankCounts(fCCContig, fCCContigID, fCCContigStart, fCCConEndPos, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        fCCContigCountOut = getFlankCounts(fCCContig, fCCContigID, fCCConStartPos, fCCContigStart, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        
        fCCResA = int2str(fCCContigCountIn) + "--" + int2str(fCCContigCountOut);

        fCCConStartPos = fCCContigEnd - fCCWinSize;
        fCCConEndPos = fCCContigEnd + fCCWinSize;
        if(fCCConStartPos < 0) fCCConStartPos = 0;
        
        fCCContigCountIn = getFlankCounts(fCCContig, fCCContigID, fCCConStartPos, fCCContigEnd, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        fCCContigCountOut = getFlankCounts(fCCContig, fCCContigID, fCCContigEnd, fCCConEndPos, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        
        fCCResA = fCCResA + "--" + int2str(fCCContigCountIn) + "--" + int2str(fCCContigCountOut);

        fCCBwaStartPos = fCCChrStart - fCCWinSize;
        fCCBwaEndPos = fCCChrStart + fCCWinSize;
        if(fCCBwaStartPos < 0) fCCBwaStartPos = 0;
        
        if(fCCOri == "forward")
        {
            fCCBwaCountIn = getFlankCounts(fCCBwa, fCCChrID, fCCChrStart, fCCBwaEndPos, fCCWinSize, (fCCBWATotal+fCCContigTotal));
            fCCBwaCountOut = getFlankCounts(fCCBwa, fCCChrID, fCCBwaStartPos, fCCChrStart, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        }
        else
        {
            fCCBwaCountIn = getFlankCounts(fCCBwa, fCCChrID, fCCBwaStartPos, fCCChrStart, fCCWinSize, (fCCBWATotal+fCCContigTotal));
            fCCBwaCountOut = getFlankCounts(fCCBwa, fCCChrID, fCCChrStart, fCCBwaEndPos, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        }
        
        fCCResB = int2str(fCCBwaCountIn) + "--" + int2str(fCCBwaCountOut);
        
        fCCBwaStartPos = fCCChrEnd - fCCWinSize;
        fCCBwaEndPos = fCCChrEnd + fCCWinSize;
        if(fCCBwaStartPos < 0) fCCBwaStartPos = 0;
        
        if(fCCOri == "forward")
        {
            fCCBwaCountIn = getFlankCounts(fCCBwa, fCCChrID, fCCBwaStartPos, fCCChrEnd, fCCWinSize, (fCCBWATotal+fCCContigTotal));
            fCCBwaCountOut = getFlankCounts(fCCBwa, fCCChrID, fCCChrEnd, fCCBwaEndPos, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        }
        else
        {
            fCCBwaCountIn = getFlankCounts(fCCBwa, fCCChrID, fCCChrEnd, fCCBwaEndPos, fCCWinSize, (fCCBWATotal+fCCContigTotal));
            fCCBwaCountOut = getFlankCounts(fCCBwa, fCCChrID, fCCBwaStartPos, fCCChrEnd, fCCWinSize, (fCCBWATotal+fCCContigTotal));
        }
        
        fCCResB = fCCResB + "--" +int2str(fCCBwaCountIn) + "--" + int2str(fCCBwaCountOut);
    }
    
    return make_pair(fCCResA, fCCResB);
}

int getBPCounts(string gBCSample, string gBCBam, string gBCID, long int gBCStart, long int gBCEnd, long int gBCPos, int gBCWin, long int gBCTotal)
{
    int gBCNum = 0;
    
    string gBCCmdTmp;
    char *gBCCmd = new char[1024];

    gBCCmdTmp = "samtools view -b -o " + gBCSample + "_GSVMiningRes/tmp.bam " + gBCBam + " " + gBCID + ":" + longint2str(gBCStart) + "-" + longint2str(gBCEnd);
    strcpy(gBCCmd, gBCCmdTmp.c_str());
    cmdExec(gBCCmd, 0);
    
    gBCCmdTmp = "samtools sort -n " + gBCSample + "_GSVMiningRes/tmp.bam -o " + gBCSample + "_GSVMiningRes/tmp.sorted.bam";
    strcpy(gBCCmd, gBCCmdTmp.c_str());
    cmdExec(gBCCmd, 0);
    
    gBCCmdTmp = "samtools view -o " + gBCSample + "_GSVMiningRes/tmp.sorted.sam " + gBCSample + "_GSVMiningRes/tmp.sorted.bam";
    strcpy(gBCCmd, gBCCmdTmp.c_str());
    cmdExec(gBCCmd, 0);
    
    string gBCFileLine;
    string gBCFileStr = gBCSample + "_GSVMiningRes/tmp.sorted.sam";
    ifstream gBCFile(gBCFileStr);
    if(!gBCFile){printf("%s \n", "Error! Cannot open tmp mapping file to read! [bpCalCov -- getBPCounts]"); exit(1);}
    
    vector<string> gBCReadA;
    vector<string> gBCReadB;
    long int gBCReadGposA = 0;
    long int gBCReadGposB = 0;
    
    gBCReadA.push_back("NA");
    gBCReadB.push_back("NA");
 
    while(1)
    {
        getline(gBCFile, gBCFileLine);
        if(gBCFile.eof()) break;
        gBCFileLine.erase(gBCFileLine.find_last_not_of(" \n\r\t")+1);
        
        gBCReadB = split(gBCFileLine, "\t");
        
        if(gBCReadA[0] == gBCReadB[0])
        {
            if(str2longint(gBCReadA[3]) < str2longint(gBCReadB[3]))
            {
                gBCReadGposA = str2longint(gBCReadA[3]);
                gBCReadGposB = str2longint(gBCReadB[3]) + gBCReadB[9].length() - 1;
            }
            else
            {
                gBCReadGposA = str2longint(gBCReadB[3]);
                gBCReadGposB = str2longint(gBCReadA[3]) + gBCReadA[9].length() - 1;
            }
            
            if((gBCReadGposA <= gBCPos)&&(gBCReadGposB >= gBCPos))
            {
                gBCNum = gBCNum + 2;
            }
        }
        else
        {
            gBCReadA = gBCReadB;
        }

    }
    
    gBCFile.close();
    
    gBCNum = gBCNum/(float(gBCWin)/1000)/(float(gBCTotal)/1000000);
    
    gBCCmdTmp = "rm " + gBCSample + "_GSVMiningRes/tmp.bam";
    strcpy(gBCCmd, gBCCmdTmp.c_str());
    cmdExec(gBCCmd, 0);
    
    gBCCmdTmp = "rm " + gBCSample + "_GSVMiningRes/tmp.sorted.bam";
    strcpy(gBCCmd, gBCCmdTmp.c_str());
    cmdExec(gBCCmd, 0);
    
    gBCCmdTmp = "rm " + gBCSample + "_GSVMiningRes/tmp.sorted.sam";
    strcpy(gBCCmd, gBCCmdTmp.c_str());
    cmdExec(gBCCmd, 0);
    
    return gBCNum;
}

pair<string, string> bpCalCov(string bCCSample, string bCCBwa, string bCCContig, int bCCGroup, string bCCContigID, string bCCChrID, long int bCCContigStart, long int bCCContigEnd, long int bCCChrStart, long int bCCChrEnd, string bCCOri, long int bCCBWATotal, long int bCCContigTotal, int bCCWinSize)
{
    long int bCCConStartPos;
    long int bCCConEndPos;
    long int bCCBwaStartPos;
    long int bCCBwaEndPos;
    
    int bCCContigCount;
    int bCCBwaCount;
    
    string bCCContigTmp;
    string bCCBwaTmp;

    if(bCCGroup == 0)
    {
        bCCConStartPos = bCCContigEnd - bCCWinSize;
        bCCConEndPos = bCCContigEnd + bCCWinSize;
        if(bCCConStartPos < 0) bCCConStartPos = 0;
        
        bCCContigCount = getBPCounts(bCCSample, bCCContig, bCCContigID, bCCConStartPos, bCCConEndPos, bCCContigEnd, bCCWinSize, (bCCBWATotal+bCCContigTotal));
        bCCContigTmp = int2str(bCCContigCount);
        
        bCCBwaStartPos = bCCChrEnd - bCCWinSize;
        bCCBwaEndPos = bCCChrEnd + bCCWinSize;
        if(bCCBwaStartPos < 0) bCCBwaStartPos = 0;
        
        bCCBwaCount = getBPCounts(bCCSample, bCCBwa, bCCChrID, bCCBwaStartPos, bCCBwaEndPos, bCCChrEnd, bCCWinSize, (bCCBWATotal+bCCContigTotal));
        bCCBwaTmp = int2str(bCCBwaCount);
    }
    else if(bCCGroup == -1)
    {
        bCCConStartPos = bCCContigStart - bCCWinSize;
        bCCConEndPos = bCCContigStart + bCCWinSize;
        if(bCCConStartPos < 0) bCCConStartPos = 0;
        
        bCCContigCount = getBPCounts(bCCSample, bCCContig, bCCContigID, bCCConStartPos, bCCConEndPos, bCCContigStart, bCCWinSize, (bCCBWATotal+bCCContigTotal));
        bCCContigTmp = int2str(bCCContigCount);
        
        bCCBwaStartPos = bCCChrStart - bCCWinSize;
        bCCBwaEndPos = bCCChrStart + bCCWinSize;
        if(bCCBwaStartPos < 0) bCCBwaStartPos = 0;
        
        bCCBwaCount = getBPCounts(bCCSample, bCCBwa, bCCChrID, bCCBwaStartPos, bCCBwaEndPos, bCCChrStart, bCCWinSize, (bCCBWATotal+bCCContigTotal));
        bCCBwaTmp = int2str(bCCBwaCount);
    }
    else
    {
        bCCConStartPos = bCCContigStart - bCCWinSize;
        bCCConEndPos = bCCContigStart + bCCWinSize;
        if(bCCConStartPos < 0) bCCConStartPos = 0;
        
        bCCContigCount = getBPCounts(bCCSample, bCCContig, bCCContigID, bCCConStartPos, bCCConEndPos, bCCContigStart, bCCWinSize, (bCCBWATotal+bCCContigTotal));
        bCCContigTmp = int2str(bCCContigCount);
        
        bCCConStartPos = bCCContigEnd - bCCWinSize;
        bCCConEndPos = bCCContigEnd + bCCWinSize;
        if(bCCConStartPos < 0) bCCConStartPos = 0;
        
        bCCContigCount = getBPCounts(bCCSample, bCCContig, bCCContigID, bCCConStartPos, bCCConEndPos, bCCContigEnd, bCCWinSize, (bCCBWATotal+bCCContigTotal));
        bCCContigTmp = bCCContigTmp + "--" + int2str(bCCContigCount);

        bCCBwaStartPos = bCCChrStart - bCCWinSize;
        bCCBwaEndPos = bCCChrStart + bCCWinSize;
        if(bCCBwaStartPos < 0) bCCBwaStartPos = 0;
        
        bCCBwaCount = getBPCounts(bCCSample, bCCBwa, bCCChrID, bCCBwaStartPos, bCCBwaEndPos, bCCChrStart, bCCWinSize, (bCCBWATotal+bCCContigTotal));
        bCCBwaTmp = int2str(bCCBwaCount);
        
        bCCBwaStartPos = bCCChrEnd - bCCWinSize;
        bCCBwaEndPos = bCCChrEnd + bCCWinSize;
        if(bCCBwaStartPos < 0) bCCBwaStartPos = 0;
        
        bCCBwaCount = getBPCounts(bCCSample, bCCBwa, bCCChrID, bCCBwaStartPos, bCCBwaEndPos, bCCChrEnd, bCCWinSize, (bCCBWATotal+bCCContigTotal));
        bCCBwaTmp = bCCBwaTmp + "--" + int2str(bCCBwaCount);
    }
    
    return make_pair(bCCContigTmp, bCCBwaTmp);
}

void flankCovFilter(string fCFID, vector<blastObj> &fCFBlast, vector<assemblyObj> &fCFAss, int fCFDist, long int fCFTotalChr, long int fCFTotalContig)
{
    cout << ">> " << currentDateTime().c_str() << " <<\t------------> Calculating Coverages of Flanking Regions\n";
    
    string bwaBam = fCFID + "_GSVMiningRes/MappingRes/" + fCFID + ".unique.sorted.bam";
    string contigBam = fCFID + "_GSVMiningRes/ReMappingRes/" + fCFID + ".remap.unique.sorted.bam";
    
    vector<blastObj> fCFBlastF;
    vector<blastObj> fCFBlastFTmp;
    vector<assemblyObj> fCFAssF;
    vector<blastObj>::iterator fCFBlastFIter;
    vector<string> fCFBlastFInfo;
    
    fCFBlast = sortBlastContig(fCFBlast);

    int fCFi = 0;
    vector<blastObj>::iterator fCFIter;
    vector<blastObj> fCFBlastList;
    regex fCFRegStr;
    
    vector<blastObj>::iterator fCFAIter;
    vector<blastObj>::iterator fCFBIter;
    
    vector<long int> fCFBlastPos;
    
    int fCFGroup = 0;
    pair<string, string> fCCCov;
    
    for(fCFi = 0; fCFi < int(fCFAss.size()); fCFi++)
    {
        fCFRegStr = fCFAss[fCFi].contigID + "\t";
        fCFBlastList = grepBlastInfo(1, fCFRegStr, fCFBlast);
        
        fCFAIter = fCFBlastList.begin();
        for(fCFIter = fCFBlastList.begin(); fCFIter != fCFBlastList.end(); fCFIter++)
        {
            fCFBIter = fCFIter;

            if(fCFAIter == fCFBlastList.begin())
            {
                fCFIter++;
                fCFBIter = fCFIter;
            }

            if((abs(fCFAIter->contigPosA - fCFBIter->contigPosA) < fCFDist)&&(abs(fCFAIter->contigPosB - fCFBIter->contigPosB) < fCFDist))
            {
                fCFBlastFTmp.push_back(*fCFAIter);
                fCFBlastPos.push_back(fCFAIter->contigPosA);
                fCFBlastPos.push_back(fCFAIter->contigPosB);
            }
            else
            {
                fCFBlastFTmp.push_back(*fCFAIter);
                fCFBlastPos.push_back(fCFAIter->contigPosA);
                fCFBlastPos.push_back(fCFAIter->contigPosB);
                
                sort(fCFBlastPos.begin(), fCFBlastPos.end());
                
                for(fCFBlastFIter = fCFBlastFTmp.begin(); fCFBlastFIter != fCFBlastFTmp.end(); fCFBlastFIter++)
                {
                    fCCCov = flankCalCov(bwaBam, contigBam, fCFGroup, fCFBlastFIter->contigID, fCFBlastFIter->chrID, fCFBlastPos.front(), fCFBlastPos.back(), fCFBlastFIter->chrPosA, fCFBlastFIter->chrPosB, fCFBlastFIter->contigOri, fCFTotalChr, fCFTotalContig, 200);
                    
                    fCFBlastFInfo = split(fCFBlastFIter->dataInfo, "\t");
                    fCFBlastFInfo[10] = fCCCov.first;
                    fCFBlastFInfo[11] = fCCCov.second;
                    fCFBlastFIter->dataInfo = joinstring(fCFBlastFInfo, "\t");
                    
                    fCFBlastF.push_back(*fCFBlastFIter);
                    
                    pair<string, string>().swap(fCCCov);
                }
                
                fCFGroup++;

                fCFBlastFTmp.clear();
                fCFBlastPos.clear();
            }
            
            fCFAIter = fCFBIter;
        }
        
        fCFBlastFTmp.push_back(*fCFAIter);
        fCFBlastPos.push_back(fCFAIter->contigPosA);
        fCFBlastPos.push_back(fCFAIter->contigPosB);
        
        sort(fCFBlastPos.begin(), fCFBlastPos.end());
        
        fCFGroup = -1;
        
        for(fCFBlastFIter = fCFBlastFTmp.begin(); fCFBlastFIter != fCFBlastFTmp.end(); fCFBlastFIter++)
        {
            fCCCov = flankCalCov(bwaBam, contigBam, fCFGroup, fCFBlastFIter->contigID, fCFBlastFIter->chrID, fCFBlastPos.front(), fCFBlastPos.back(), fCFBlastFIter->chrPosA, fCFBlastFIter->chrPosB, fCFBlastFIter->contigOri, fCFTotalChr, fCFTotalContig, 200);
            
            fCFBlastFInfo = split(fCFBlastFIter->dataInfo, "\t");
            fCFBlastFInfo[10] = fCCCov.first;
            fCFBlastFInfo[11] = fCCCov.second;
            fCFBlastFIter->dataInfo = joinstring(fCFBlastFInfo, "\t");
            
            fCFBlastF.push_back(*fCFBlastFIter);
            
            pair<string, string>().swap(fCCCov);
        }
        
        fCFBlastFTmp.clear();
        fCFBlastPos.clear();
        fCFBlastList.clear();
        
        fCFGroup = 0;
    }
    
    cout << ">> " << currentDateTime().c_str() << " <<\t------------> Filtering 0 Data\n";
    vector<string> fCFBlastFCovInfo;
    
    fCFBlastFTmp = fCFBlastF;
    fCFBlastF.clear();
    for(fCFBlastFIter = fCFBlastFTmp.begin(); fCFBlastFIter != fCFBlastFTmp.end(); fCFBlastFIter++)
    {
        fCFBlastFInfo = split(fCFBlastFIter->dataInfo, "\t");
        fCFBlastFCovInfo = split(fCFBlastFInfo[10], "--");
        
        if(fCFBlastFCovInfo.size() == 2)
        {
            if((str2int(fCFBlastFCovInfo[0]) != 0)&&(str2int(fCFBlastFCovInfo[1]) != 0)){fCFBlastF.push_back(*fCFBlastFIter);}
        }
        else
        {
            if(((str2int(fCFBlastFCovInfo[0]) != 0)&&(str2int(fCFBlastFCovInfo[1]) != 0))||((str2int(fCFBlastFCovInfo[2]) != 0)&&(str2int(fCFBlastFCovInfo[3]) != 0))){fCFBlastF.push_back(*fCFBlastFIter);}
        }
    }
    
    fCFBlastFTmp.clear();
    fCFBlastFInfo.clear();
    fCFBlastFCovInfo.clear();
    
    cout << ">> " << currentDateTime().c_str() << " <<\t------------> Filtering Low Data\n";
    
    fCFBlastFTmp = fCFBlastF;
    fCFBlastF.clear();
    
    vector<string> fCFBlastFCovBwaInfo;
    
    float covCmpA;
    float covCmpB;
    float covBwaCmpA;
    float covBwaCmpB;
    
    for(fCFBlastFIter = fCFBlastFTmp.begin(); fCFBlastFIter != fCFBlastFTmp.end(); fCFBlastFIter++)
    {
        fCFBlastFInfo = split(fCFBlastFIter->dataInfo, "\t");
        fCFBlastFCovInfo = split(fCFBlastFInfo[10], "--");
        fCFBlastFCovBwaInfo = split(fCFBlastFInfo[11], "--");
        
        if(fCFBlastFCovInfo.size() == 2)
        {
            covCmpA = (str2int(fCFBlastFCovInfo[0]) + str2int(fCFBlastFCovInfo[1]))/2;
            covBwaCmpA = (str2int(fCFBlastFCovBwaInfo[0]) + str2int(fCFBlastFCovBwaInfo[1]))/2;
            
            if((covBwaCmpA/covCmpA) < 20){fCFBlastF.push_back(*fCFBlastFIter);}
        }
        else
        {
            covCmpA = (str2int(fCFBlastFCovInfo[0]) + str2int(fCFBlastFCovInfo[1]))/2;
            covBwaCmpA = (str2int(fCFBlastFCovBwaInfo[0]) + str2int(fCFBlastFCovBwaInfo[1]))/2;
            
            covCmpB = (str2int(fCFBlastFCovInfo[2]) + str2int(fCFBlastFCovInfo[3]))/2;
            covBwaCmpB = (str2int(fCFBlastFCovBwaInfo[2]) + str2int(fCFBlastFCovBwaInfo[3]))/2;
            
            if((covCmpA != 0)&&(covCmpB != 0))
            {
                if(((covBwaCmpA/covCmpA) < 20)||((covBwaCmpB/covCmpB) < 20)){fCFBlastF.push_back(*fCFBlastFIter);}
            }
            else if(covCmpA == 0)
            {
                if((covBwaCmpB/covCmpB) < 20){fCFBlastF.push_back(*fCFBlastFIter);}
            }
            else if(covCmpB == 0)
            {
                if((covBwaCmpA/covCmpA) < 20){fCFBlastF.push_back(*fCFBlastFIter);}
            }
        }
    }
    
    fCFBlastFTmp.clear();
    fCFBlastFInfo.clear();
    fCFBlastFCovInfo.clear();
    fCFBlastFCovBwaInfo.clear();
    
    cout << ">> " << currentDateTime().c_str() << " <<\t------------> Filtering Inconsecutive Data (Incomplete Alignment)\n";
    
    fCFGroup = 0;
    
    fCFBlastFTmp = fCFBlastF;
    fCFBlastF.clear();
    fCFBlastList.clear();
    
    vector<blastObj> fCFFragList;
    vector<string> blastContigID;
    
    for(fCFBlastFIter = fCFBlastFTmp.begin(); fCFBlastFIter != fCFBlastFTmp.end(); fCFBlastFIter++)
    {
        blastContigID.push_back(fCFBlastFIter->contigID);
    }
    
    blastContigID.erase(unique(blastContigID.begin(), blastContigID.end()), blastContigID.end());

    for(fCFi = 0; fCFi < int(blastContigID.size()); fCFi++)
    {
        fCFRegStr = blastContigID[fCFi] + "\t";
        fCFBlastList = grepBlastInfo(1, fCFRegStr, fCFBlastFTmp);
        
        if(fCFBlastList.size() > 1)
        {
            fCFAIter = fCFBlastList.begin();
            for(fCFIter = fCFBlastList.begin(); fCFIter != fCFBlastList.end(); fCFIter++)
            {
                fCFBIter = fCFIter;
                
                if(fCFAIter == fCFBlastList.begin())
                {
                    fCFIter++;
                    fCFBIter = fCFIter;
                }
                
                if((abs(fCFAIter->contigPosA - fCFBIter->contigPosA) < fCFDist)||(abs(fCFAIter->contigPosB - fCFBIter->contigPosB) < fCFDist))
                {
                    fCFFragList.push_back(*fCFAIter);
                }
                else
                {
                    fCFFragList.push_back(*fCFAIter);
                    
                    fCFGroup++;
                    fCFBlastF.insert(fCFBlastF.end(), fCFFragList.begin(), fCFFragList.end());
                    
                    fCFFragList.clear();
                }
                
                fCFAIter = fCFBIter;
            }
            
            fCFFragList.push_back(*fCFAIter);
            if(fCFGroup != 0){fCFBlastF.insert(fCFBlastF.end(), fCFFragList.begin(), fCFFragList.end());}
            
            fCFFragList.clear();
            fCFBlastList.clear();
            fCFGroup = 0;
        }
    }
    
    cout << ">> " << currentDateTime().c_str() << " <<\t------------> Generating Positions\n";
    
    string fCFFileStr = fCFID + "_GSVMiningRes/BlastRes/" + fCFID + ".GRearr.candidates.bpPos.txt";
    ofstream fCFBlastOut(fCFFileStr);
    if(!fCFBlastOut){printf("%s \n", "Error! Cannot open file to write! [flankCovFilter -- Generating Positions]"); exit(1);}
    
    fCFGroup = 0;
    
    fCFBlastFTmp = fCFBlastF;
    fCFBlastList.clear();
    fCFFragList.clear();
    blastContigID.clear();
    
    for(fCFBlastFIter = fCFBlastFTmp.begin(); fCFBlastFIter != fCFBlastFTmp.end(); fCFBlastFIter++)
    {
        blastContigID.push_back(fCFBlastFIter->contigID);
    }
    
    blastContigID.erase(unique(blastContigID.begin(), blastContigID.end()), blastContigID.end());
    
    for(fCFi = 0; fCFi < int(blastContigID.size()); fCFi++)
    {
        fCFRegStr = blastContigID[fCFi] + "\t";
        fCFBlastList = grepBlastInfo(1, fCFRegStr, fCFBlastFTmp);
        
        fCFAIter = fCFBlastList.begin();
        for(fCFIter = fCFBlastList.begin(); fCFIter != fCFBlastList.end(); fCFIter++)
        {
            fCFBIter = fCFIter;
            
            if(fCFAIter == fCFBlastList.begin())
            {
                fCFIter++;
                fCFBIter = fCFIter;
            }
            
            if((abs(fCFAIter->contigPosA - fCFBIter->contigPosA) < fCFDist)||(abs(fCFAIter->contigPosB - fCFBIter->contigPosB) < fCFDist))
            {
                fCFFragList.push_back(*fCFAIter);
            }
            else
            {
                fCFFragList.push_back(*fCFAIter);
                
                for(fCFBlastFIter = fCFFragList.begin(); fCFBlastFIter != fCFFragList.end(); fCFBlastFIter++)
                {
                    if(fCFGroup == 0)
                    {
                        fCFBlastOut << fCFBlastFIter->contigID << "\t" << fCFBlastFIter->chrID << "\t" << fCFBlastFIter->chrPosB << "\n";
                    }
                    else
                    {
                        fCFBlastOut << fCFBlastFIter->contigID << "\t" << fCFBlastFIter->chrID << "\t" << fCFBlastFIter->chrPosA << "\n";
                        fCFBlastOut << fCFBlastFIter->contigID << "\t" << fCFBlastFIter->chrID << "\t" << fCFBlastFIter->chrPosB << "\n";
                    }
                }
                
                fCFGroup++;
                
                fCFFragList.clear();
            }
            
            fCFAIter = fCFBIter;
        }
        
        fCFFragList.push_back(*fCFAIter);
        
        for(fCFBlastFIter = fCFFragList.begin(); fCFBlastFIter != fCFFragList.end(); fCFBlastFIter++)
        {
            fCFBlastOut << fCFBlastFIter->contigID << "\t" << fCFBlastFIter->chrID << "\t" << fCFBlastFIter->chrPosA << "\n";
        }
        
        fCFFragList.clear();
        fCFBlastList.clear();
        fCFGroup = 0;
    }
    
    fCFBlastOut.close();
    
    fCFBlast.clear();
    fCFBlast = fCFBlastF;
    
    vector<blastObj>().swap(fCFBlastF);
    vector<blastObj>().swap(fCFBlastFTmp);
    vector<blastObj>().swap(fCFBlastList);
    vector<assemblyObj>().swap(fCFAssF);
}

void bpCovFilter(string bCFID, vector<blastObj> &bCFBlast, vector<assemblyObj> &bCFAss, int bCFDist, int bCFMaxLen, long int bCFTotalChr, long int bCFTotalContig)
{
    cout << ">> " << currentDateTime().c_str() << " <<\t------------> Calculating Coverages of Breakpoint Regions\n";
    
    string bwaBam = bCFID + "_GSVMiningRes/MappingRes/" + bCFID + ".unique.sorted.bam";
    string contigBam = bCFID + "_GSVMiningRes/ReMappingRes/" + bCFID + ".remap.unique.sorted.bam";
    
    vector<blastObj> bCFBlastF;
    vector<blastObj> bCFBlastFTmp;
    vector<assemblyObj> bCFAssF;
    vector<blastObj>::iterator bCFBlastFIter;
    vector<string> bCFBlastFInfo;
    vector<string> blastContigID;
    
    bCFBlast = sortBlastContig(bCFBlast);
    
    int bCFi = 0;
    vector<blastObj>::iterator bCFIter;
    vector<blastObj> bCFBlastList;
    regex bCFRegStr;
    
    vector<blastObj>::iterator bCFAIter;
    vector<blastObj>::iterator bCFBIter;
    
    vector<long int> bCFBlastPos;
    
    int bCFGroup = 0;
    pair<string, string> bCCCov;
    
    for(bCFBlastFIter = bCFBlast.begin(); bCFBlastFIter != bCFBlast.end(); bCFBlastFIter++)
    {
        blastContigID.push_back(bCFBlastFIter->contigID);
    }
    
    blastContigID.erase(unique(blastContigID.begin(), blastContigID.end()), blastContigID.end());
    
    for(bCFi = 0; bCFi < int(blastContigID.size()); bCFi++)
    {
        bCFRegStr = blastContigID[bCFi] + "\t";
        bCFBlastList = grepBlastInfo(1, bCFRegStr, bCFBlast);
        
        bCFAIter = bCFBlastList.begin();
        for(bCFIter = bCFBlastList.begin(); bCFIter != bCFBlastList.end(); bCFIter++)
        {
            bCFBIter = bCFIter;
            
            if(bCFAIter == bCFBlastList.begin())
            {
                bCFIter++;
                bCFBIter = bCFIter;
            }
            
            if((abs(bCFAIter->contigPosA - bCFBIter->contigPosA) < bCFDist)&&(abs(bCFAIter->contigPosB - bCFBIter->contigPosB) < bCFDist))
            {
                bCFBlastFTmp.push_back(*bCFAIter);
                bCFBlastPos.push_back(bCFAIter->contigPosA);
                bCFBlastPos.push_back(bCFAIter->contigPosB);
            }
            else
            {
                bCFBlastFTmp.push_back(*bCFAIter);
                bCFBlastPos.push_back(bCFAIter->contigPosA);
                bCFBlastPos.push_back(bCFAIter->contigPosB);
                
                sort(bCFBlastPos.begin(), bCFBlastPos.end());
                
                for(bCFBlastFIter = bCFBlastFTmp.begin(); bCFBlastFIter != bCFBlastFTmp.end(); bCFBlastFIter++)
                {
                    bCCCov = bpCalCov(bCFID, bwaBam, contigBam, bCFGroup, bCFBlastFIter->contigID, bCFBlastFIter->chrID, bCFBlastPos.front(), bCFBlastPos.back(), bCFBlastFIter->chrPosA, bCFBlastFIter->chrPosB, bCFBlastFIter->contigOri, bCFTotalChr, bCFTotalContig, bCFMaxLen);
                    
                    bCFBlastFInfo = split(bCFBlastFIter->dataInfo, "\t");
                    bCFBlastFInfo[4] = bCCCov.first;
                    bCFBlastFInfo[5] = bCCCov.second;
                    bCFBlastFIter->dataInfo = joinstring(bCFBlastFInfo, "\t");
                    
                    bCFBlastF.push_back(*bCFBlastFIter);
                    
                    pair<string, string>().swap(bCCCov);
                }
                
                bCFGroup++;
                
                bCFBlastFTmp.clear();
                bCFBlastPos.clear();
            }
            
            bCFAIter = bCFBIter;
        }
        
        bCFBlastFTmp.push_back(*bCFAIter);
        bCFBlastPos.push_back(bCFAIter->contigPosA);
        bCFBlastPos.push_back(bCFAIter->contigPosB);
        
        sort(bCFBlastPos.begin(), bCFBlastPos.end());
        
        bCFGroup = -1;
        
        for(bCFBlastFIter = bCFBlastFTmp.begin(); bCFBlastFIter != bCFBlastFTmp.end(); bCFBlastFIter++)
        {
            bCCCov = bpCalCov(bCFID, bwaBam, contigBam, bCFGroup, bCFBlastFIter->contigID, bCFBlastFIter->chrID, bCFBlastPos.front(), bCFBlastPos.back(), bCFBlastFIter->chrPosA, bCFBlastFIter->chrPosB, bCFBlastFIter->contigOri, bCFTotalChr, bCFTotalContig, bCFMaxLen);
            
            bCFBlastFInfo = split(bCFBlastFIter->dataInfo, "\t");
            bCFBlastFInfo[4] = bCCCov.first;
            bCFBlastFInfo[5] = bCCCov.second;
            bCFBlastFIter->dataInfo = joinstring(bCFBlastFInfo, "\t");
            
            bCFBlastF.push_back(*bCFBlastFIter);
            
            pair<string, string>().swap(bCCCov);
        }
        
        bCFBlastFTmp.clear();
        bCFBlastPos.clear();
        bCFBlastList.clear();
        
        bCFGroup = 0;
    }
    
    bCFBlast.clear();
    bCFBlast = bCFBlastF;
    
    vector<blastObj>().swap(bCFBlastF);
    vector<blastObj>().swap(bCFBlastFTmp);
    vector<blastObj>().swap(bCFBlastList);
    vector<assemblyObj>().swap(bCFAssF);
}

void GSVCovFilter(string GSVCFID, int GSVCFDist, int GSVCFLen, long int GSVCFTotalChr, long int GSVCFTotalContig)
{
    string GSVCFLine;
    string GSVCFFileStr;
    vector<string> GSVCFDataStr;
    
    blastObj GSVCFBlastObj;
    vector<blastObj> GSVCFBlastRes;
    GSVCFFileStr = GSVCFID + "_GSVMiningRes/BlastRes/" + GSVCFID + ".GRearr.candidates.filter.rm.txt";
    ifstream GSVCFFileBlast(GSVCFFileStr);
    if(!GSVCFFileBlast){printf("%s \n", "Error! Cannot open blast result list file to read! [GSVPrimaryFilter]"); exit(1);}
    while(1)
    {
        getline(GSVCFFileBlast, GSVCFLine); // careful! EOF check must be after this cmd as otherwise cannot get right EOF position.
        if(GSVCFFileBlast.eof()) break;
        GSVCFLine.erase(GSVCFLine.find_last_not_of(" \n\r\t")+1);
        
        GSVCFDataStr = split(GSVCFLine, "\t");
        GSVCFBlastObj.getBlastInfo(GSVCFDataStr[0], GSVCFDataStr[1], str2int(GSVCFDataStr[3]), str2int(GSVCFDataStr[6]), str2int(GSVCFDataStr[7]), str2int(GSVCFDataStr[8]), str2int(GSVCFDataStr[9]), GSVCFLine);
        GSVCFBlastObj.contigOri = GSVCFDataStr[12];
        GSVCFBlastRes.push_back(GSVCFBlastObj);
    }
    
    GSVCFFileBlast.close();
    
    assemblyObj GSVCFAssObj;
    vector<assemblyObj> GSVCFAssRes;
    GSVCFFileStr = GSVCFID + "_GSVMiningRes/BlastRes/" + GSVCFID + ".GRearr.candidates.filter.rm.rpkm.txt";
    ifstream GSVCFFileAss(GSVCFFileStr);
    if(!GSVCFFileAss){printf("%s \n", "Error! Cannot open blast result rpkm file to read! [GSVPrimaryFilter]"); exit(1);}
    
    while(1)
    {
        getline(GSVCFFileAss, GSVCFLine); // careful! EOF check must be after this cmd as otherwise cannot get right EOF position.
        if(GSVCFFileAss.eof()) break;
        GSVCFLine.erase(GSVCFLine.find_last_not_of(" \n\r\t")+1);
        
        GSVCFDataStr = split(GSVCFLine, "\t");
        GSVCFAssObj.getAssInfo(GSVCFDataStr[0], str2int(GSVCFDataStr[1]), str2longint(GSVCFDataStr[2]), str2float(GSVCFDataStr[3]), str2float(GSVCFDataStr[4]));
        GSVCFAssRes.push_back(GSVCFAssObj);
    }
    
    GSVCFFileAss.close();

    cout << ">> " << currentDateTime().c_str() << " <<\t--------> Flanking Coverage Filtering\n";
    flankCovFilter(GSVCFID, GSVCFBlastRes, GSVCFAssRes, GSVCFDist, GSVCFTotalChr, GSVCFTotalContig);
    
    cout << ">> " << currentDateTime().c_str() << " <<\t--------> Breakpoint Coverage Filtering\n";
    bpCovFilter(GSVCFID, GSVCFBlastRes, GSVCFAssRes, GSVCFDist, 1200, GSVCFTotalChr, GSVCFTotalContig);
    
    GSVCFFileStr = GSVCFID + "_GSVMiningRes/BlastRes/" + GSVCFID + ".GRearr.candidates.bpCov.txt";
    ofstream GSVCFFileBlastOut(GSVCFFileStr);
    if(!GSVCFFileBlastOut){printf("%s \n", "Error! Cannot open file to write! [flankCovFilter -- Generating Positions]"); exit(1);}
    
    vector<blastObj>::iterator GSVCFIter;
    
    for(GSVCFIter = GSVCFBlastRes.begin(); GSVCFIter != GSVCFBlastRes.end(); GSVCFIter++)
    {
        GSVCFFileBlastOut << GSVCFIter->dataInfo << "\n";
    }
    
    GSVCFFileBlastOut.close();
    
    // free vector
    vector<string>().swap(GSVCFDataStr);

}
