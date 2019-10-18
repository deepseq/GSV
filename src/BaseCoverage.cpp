/*
 
 BaseCoverage.cpp
 
 This program is to do:
 1. Count base coverage for each chromosome
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
*/

#include "BaseCoverage.h"
#include "CurrentTime.h"
#include "CmdExecution.h"

#include <cstring>
#include <iostream>
#include <vector>

using namespace std;

void baseCov(string bCID, vector<string> bCChr)
{
    vector<string>::iterator bCIter;
    
    string baseCovTmp;
    char *baseCovCmd = new char[256];
    
    for(bCIter = bCChr.begin(); bCIter != bCChr.end(); bCIter++)
    {
        cout << ">> " << currentDateTime().c_str() << " <<\t--------> " << *bCIter << "\n";
        
        baseCovTmp = "samtools view -b -o " + bCID + "_GSVMiningRes/MappingRes/" + bCID + "." + *bCIter + ".bam " + bCID + "_GSVMiningRes/MappingRes/" + bCID + ".unique.sorted.bam " + *bCIter;
        strcpy(baseCovCmd, baseCovTmp.c_str());
        cmdExec(baseCovCmd, 0);

        baseCovTmp = "bedtools genomecov -ibam " + bCID + "_GSVMiningRes/MappingRes/" + bCID + "." + *bCIter + ".bam -d > " + bCID + "_GSVMiningRes/MappingRes/" + bCID + "." + *bCIter + ".basecov.txt";
        strcpy(baseCovCmd, baseCovTmp.c_str());
        cmdExec(baseCovCmd, 0);
        
        baseCovTmp = "rm " + bCID + "_GSVMiningRes/MappingRes/" + bCID + "." + *bCIter + ".bam";
        strcpy(baseCovCmd, baseCovTmp.c_str());
        cmdExec(baseCovCmd, 0);
    }
    
    delete[] baseCovCmd;
}
