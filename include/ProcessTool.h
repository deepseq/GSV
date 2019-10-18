//
//  ProcessTool.h
//  GSVMining
//
//  Created by admin on 14-12-18.
//  Copyright (c) 2014年 University of Nottingham. All rights reserved.
//

#ifndef __GSVMining__ProcessTool__
#define __GSVMining__ProcessTool__

#include <iostream>
#include <regex>
#include <vector>
#include <algorithm>

using namespace std;

class blastObj
{
public:
    string contigID;
    string chrID;
    int contigAlign;
    int contigPosA;
    int contigPosB;
    int chrPosA;
    int chrPosB;
    string dataInfo;
    
    string contigOri;
    
public:
    void getBlastInfo(string subConID, string subChrID, int subConAlign, int subConPosA, int subConPosB, int subChrPosA, int subChrPosB, string subData);
};

class assemblyObj
{
public:
    string contigID;
    string contigSeq;
    int contigLen;
    long int contigCount;
    float contigRPKM;
    float contigLogRPKM;
    
public:
    void getAssInfo(string subConID, int subConLen, long int subConCount, float subConRPKM, float subConLogRPKM);
};

vector<blastObj> grepBlastInfo(int gBIType, regex gBIKey, vector<blastObj> gBIIn);
vector<blastObj> sortBlastContig(vector<blastObj> sBConIn);
vector<blastObj> sortBlastChr(vector<blastObj> sBChrIn);
vector<blastObj> uniqueBlast(vector<blastObj> uBIn);

#endif /* defined(__GSVMining__ProcessTool__) */
