/*
 
 DataPreparation.cpp
  
 This program is to do:
 1. create folders
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
*/

#include "DataPreparation.h"
#include "CmdExecution.h"

#include <cstring>
#include <iostream>
#include <sys/stat.h>

using namespace std;

void dataPre(string dPID, string dPRef)
{
    string dataPreTmp;
    char *dataPreCmd = new char[1024];
    struct stat dataPreSb;
    
    dataPreTmp = dPID + "_GSVMiningRes";
    strcpy(dataPreCmd, dataPreTmp.c_str());
    
    // check if tmp directory exists
    if(stat(dataPreCmd, &dataPreSb) == 0 && S_ISDIR(dataPreSb.st_mode))
    {
        dataPreTmp = "rm -rf " + dPID + "_GSVMiningRes";
        strcpy(dataPreCmd, dataPreTmp.c_str());
        cmdExec(dataPreCmd, 0);
    }
    
    dataPreTmp = "mkdir " + dPID + "_GSVMiningRes";
    strcpy(dataPreCmd, dataPreTmp.c_str());
    cmdExec(dataPreCmd, 0);
    
    dataPreTmp = "mkdir " + dPID + "_GSVMiningRes/MappingRes";
    strcpy(dataPreCmd, dataPreTmp.c_str());
    cmdExec(dataPreCmd, 0);
    
    dataPreTmp = "mkdir " + dPID + "_GSVMiningRes/AssemblyRes";
    strcpy(dataPreCmd, dataPreTmp.c_str());
    cmdExec(dataPreCmd, 0);

    dataPreTmp = "mkdir " + dPID + "_GSVMiningRes/BlastRes";
    strcpy(dataPreCmd, dataPreTmp.c_str());
    cmdExec(dataPreCmd, 0);
    
    dataPreTmp = "mkdir " + dPID + "_GSVMiningRes/ReMappingRes";
    strcpy(dataPreCmd, dataPreTmp.c_str());
    cmdExec(dataPreCmd, 0);

    dataPreTmp = "mkdir " + dPID + "_GSVMiningRes/UnmappedReads";
    strcpy(dataPreCmd, dataPreTmp.c_str());
    cmdExec(dataPreCmd, 0);
    
    delete[] dataPreCmd;
}

