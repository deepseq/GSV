/*
 
 MappingSummary.cpp
 
 This program is to do:
 1. get unmapped reads
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
*/

#include "CmdExecution.h"
#include "CurrentTime.h"
#include "MappingSummary.h"
#include "StringOperation.h"

#include <iostream>
#include <fstream>
#include <regex>
#include <vector>

using namespace std;

///////////////////
/* sub functions */
///////////////////

string DecToBin(int DBNum)
{
    if(DBNum == 0){return "0";}
    if(DBNum == 1){return "1";}
    
    if(DBNum%2 == 0)
    {
        return DecToBin(DBNum/2) + "0";
    }
    else
    {
        return DecToBin(DBNum/2) + "1";
    }
}

int BinToDec(string BDNum)
{
    int BDRes = 0;
    int BDPow = 1;
    int i = 0;
    
    for(i = int(BDNum.length()-1); i >= 0; --i, BDPow <<= 1)
    {
        BDRes += (BDNum[i] - '0') * BDPow;
    }
    
    return BDRes;
}

////////////////////
/* main functions */
////////////////////

// 147v99 or 83v163 mapped with right orientation and right insert size
// 73v133 or 89v121 or 165v181 or 101v117 or 153v185 or 69v137 one of mate is unmapped
// 77v141 both unmapped
// 67v131 or 115v179 mapped with wrong orientation
// 81v161 or 97v145 or 65v129 or 113v177 mapped with wrong insert size, so possibly multiple mapping locus

long int mapSum(string mSID)
{
    string mapSumLine;
    string mapSumFH;
    
    string mapSumTmp;
    char *mapSumCmd = new char[1024];
    
    regex mSRegStr("@[A-Z][A-Z]\t");
    
    // claim variables
    string mSTotalMappedStr;
    long int mSTotalMapped = 0;
    
    /* get unmapped reads */
    
    mapSumTmp = "samtools view -H " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".bam > " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".header.sam";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "samtools view -bh -F 2 -o " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".unmapped.bam " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".bam";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "samtools view -c " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".unique.sorted.bam";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    mSTotalMappedStr = cmdGetOut(mapSumCmd, 0);
    mSTotalMapped = str2int(mSTotalMappedStr);
    
    mapSumTmp = "bam2fastx -q -A -o " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + " -P -N " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".unmapped.bam";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "mv " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".1. " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".unmapped.tmp.read_1.fastq";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "mv " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".2. " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".unmapped.tmp.read_2.fastq";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    
    /* get clipped reads */
    
    mapSumTmp = "samtools view " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".unique.bam | awk '{if(($6~/S/)&&($0~/XP:Z/)){print $1}}' > " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".XP.readID.txt | sort | uniq";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "awk 'NR==FNR{a[$1]=$1;next}{if(a[$1]){print $0}}' " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".XP.readID.txt " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".unique.sam > " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".XP.paired.tmp.sam";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "cat " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".header.sam " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".XP.paired.tmp.sam > " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".XP.paired.sam";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "bam2fastx -q -s -A -o " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + " -P -N " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".XP.paired.sam";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "mv " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".1. " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".clipped.read_1.fastq";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "mv " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".2. " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".clipped.read_2.fastq";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    /* merge reads */
    
    mapSumTmp = "cat " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".unmapped.tmp.read_1.fastq " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".clipped.read_1.fastq > " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".unmapped.read_1.fastq";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "cat " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".unmapped.tmp.read_2.fastq " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".clipped.read_2.fastq > " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".unmapped.read_2.fastq";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    /* rm files */
    
    //move to MappingReads.cpp
    //mapSumTmp = "rm -rf " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".sam";
    //strcpy(mapSumCmd, mapSumTmp.c_str());
    //cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "rm -rf " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".unique.bam";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "rm -rf " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".unique.sam";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "rm -rf " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".unmapped.sam";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "rm -rf " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".XP.readID.txt";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "rm -rf " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".XP.paired.tmp.sam";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "rm -rf " + mSID + "_GSVMiningRes/MappingRes/" + mSID + ".XP.paired.sam";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "rm -rf " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".unmapped.tmp.read_*.fastq";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    mapSumTmp = "rm -rf " + mSID + "_GSVMiningRes/UnmappedReads/" + mSID + ".clipped.read_*.fastq";
    strcpy(mapSumCmd, mapSumTmp.c_str());
    cmdExec(mapSumCmd, 0);
    
    // free memory
    
    delete[] mapSumCmd;
    
    return mSTotalMapped;
}