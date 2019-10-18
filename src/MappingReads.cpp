/*
 
 MappingReads.cpp
 
 This program is to do:
 1. map reads to reference
 2. filter mapping results
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
*/

#include "MappingReads.h"
#include "CmdExecution.h"
#include "StringOperation.h"

#include <cstring>
#include <iostream>

using namespace std;

void mapReads(string mRRef, string mRDir, string mRID, int mRThread)
{
    string mapReadsLogStr = mRID + ".log.txt";
    
    string mapReadsTmp;
    char *mapReadsCmd = new char[1024];
    
    mapReadsTmp = "bwa mem -t " + int2str(mRThread) + " " + mRRef + " " + mRDir + "/" + mRID + ".read_1.fastq.gz " + mRDir + "/" + mRID + ".read_2.fastq.gz > " + mRID + "_GSVMiningRes/MappingRes/" + mRID + ".sam 2>> " + mapReadsLogStr;
    strcpy(mapReadsCmd, mapReadsTmp.c_str());
    cmdExec(mapReadsCmd, 0);
    
    mapReadsTmp = "samtools view -bS " + mRID + "_GSVMiningRes/MappingRes/" + mRID + ".sam > " + mRID + "_GSVMiningRes/MappingRes/" + mRID + ".bam";
    strcpy(mapReadsCmd, mapReadsTmp.c_str());
    cmdExec(mapReadsCmd, 0);
    
    mapReadsTmp = "samtools flagstat " + mRID + "_GSVMiningRes/MappingRes/" + mRID + ".bam > " + mRID + "_GSVMiningRes/MappingRes/" + mRID + ".stat.txt";
    strcpy(mapReadsCmd, mapReadsTmp.c_str());
    cmdExec(mapReadsCmd, 0);
    
    mapReadsTmp = "awk '{if(($2==99)||($2==147)||($2==83)||($2==163)||($1~/^@/)){print $0}}' " + mRID + "_GSVMiningRes/MappingRes/" + mRID + ".sam | grep \"XS:i:0\\|^@\" > " + mRID + "_GSVMiningRes/MappingRes/" + mRID + ".unique.sam";
    strcpy(mapReadsCmd, mapReadsTmp.c_str());
    cmdExec(mapReadsCmd, 0);
    
    mapReadsTmp = "samtools view -bS " + mRID + "_GSVMiningRes/MappingRes/" + mRID + ".unique.sam > " + mRID + "_GSVMiningRes/MappingRes/" + mRID + ".unique.bam";
    strcpy(mapReadsCmd, mapReadsTmp.c_str());
    cmdExec(mapReadsCmd, 0);

    mapReadsTmp = "samtools sort -@ " + int2str(mRThread) + " " + mRID + "_GSVMiningRes/MappingRes/" + mRID + ".unique.bam -o " + mRID + "_GSVMiningRes/MappingRes/" + mRID + ".unique.sorted.bam";
    strcpy(mapReadsCmd, mapReadsTmp.c_str());
    cmdExec(mapReadsCmd, 0);
    
    mapReadsTmp = "samtools index " + mRID + "_GSVMiningRes/MappingRes/" + mRID + ".unique.sorted.bam";
    strcpy(mapReadsCmd, mapReadsTmp.c_str());
    cmdExec(mapReadsCmd, 0);

    mapReadsTmp = "samtools flagstat " + mRID + "_GSVMiningRes/MappingRes/" + mRID + ".unique.sorted.bam > " + mRID + "_GSVMiningRes/MappingRes/" + mRID + ".unique.stat.txt";
    strcpy(mapReadsCmd, mapReadsTmp.c_str());
    cmdExec(mapReadsCmd, 0);

    /* rm files */
    
    mapReadsTmp = "rm -rf " + mRID + "_GSVMiningRes/MappingRes/" + mRID + ".sam";
    strcpy(mapReadsCmd, mapReadsTmp.c_str());
    cmdExec(mapReadsCmd, 0);
    
    delete[] mapReadsCmd;
}
