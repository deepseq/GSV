//
//  ProcessTool.cpp
//  GSVMining
//
//  Created by admin on 14-12-18.
//  Copyright (c) 2014年 University of Nottingham. All rights reserved.
//

#include "ProcessTool.h"
#include "StringOperation.h"

#include <regex>
#include <vector>
#include <algorithm>

using namespace std;

///////////////////
/* program class */
///////////////////

void blastObj::getBlastInfo(string subConID, string subChrID, int subConAlign, int subConPosA, int subConPosB, int subChrPosA, int subChrPosB, string subData)
{
    contigID = subConID;
    chrID = subChrID;
    contigAlign = subConAlign;
    contigPosA = subConPosA;
    contigPosB = subConPosB;
    chrPosA = subChrPosA;
    chrPosB = subChrPosB;
    dataInfo = subData;
}

void assemblyObj::getAssInfo(string subConID, int subConLen, long subConCount, float subConRPKM, float subConLogRPKM)
{
    contigID = subConID;
    contigLen = subConLen;
    contigCount = subConCount;
    contigRPKM = subConRPKM;
    contigLogRPKM = subConLogRPKM;
}

////////////////////
/* tool functions */
////////////////////

vector<blastObj> grepBlastInfo(int gBIType, regex gBIKey, vector<blastObj> gBIIn)
{
    vector<blastObj> gBIRes;
    int gBIi;
    
    for(gBIi = 0; gBIi < int(gBIIn.size()); gBIi++)
    {
        if(gBIType == 1)
        {
            if(regex_search(gBIIn[gBIi].dataInfo, gBIKey, regex_constants::match_continuous))
            {
                gBIRes.push_back(gBIIn[gBIi]);
            }
        }
        else
        {
            if(regex_search(gBIIn[gBIi].dataInfo, gBIKey))
            {
                gBIRes.push_back(gBIIn[gBIi]);
            }
        }
    }
    
    return gBIRes;
}

bool cmpBlastContig(blastObj cBConA, blastObj cBConB)
{
    if(cBConA.contigID.compare(cBConB.contigID) < 0) return true;
    if(cBConA.contigID.compare(cBConB.contigID) == 0) return (cBConA.contigPosA < cBConB.contigPosA);
    
    return false;
}

vector<blastObj> sortBlastContig(vector<blastObj> sBConIn)
{
    sort(sBConIn.begin(), sBConIn.end(), cmpBlastContig);
    return sBConIn;
}

bool cmpBlastChr(blastObj cBChrA, blastObj cBChrB)
{
    if(cBChrA.chrID.compare(cBChrB.chrID) < 0) return true;
    if(cBChrA.chrID.compare(cBChrB.chrID) == 0) return (cBChrA.chrPosA < cBChrB.chrPosA);
    
    return false;
}

vector<blastObj> sortBlastChr(vector<blastObj> sBChrIn)
{
    sort(sBChrIn.begin(), sBChrIn.end(), cmpBlastChr);
    return sBChrIn;
}

vector<blastObj> uniqueBlast(vector<blastObj> uBIn)
{
    vector<blastObj> uBOut;
    
    int uBi;
    vector<string> uBInfo;
    for(uBi = 0; uBi < int(uBIn.size()); uBi++)
    {
        uBInfo.push_back(uBIn[uBi].dataInfo);
    }
    
    uBInfo.erase(unique(uBInfo.begin(), uBInfo.end()), uBInfo.end());
    
    vector<string> uBInfoStr;
    blastObj uBOutTmp;
    for(uBi = 0; uBi < int(uBInfo.size()); uBi++)
    {
        uBInfoStr = split(uBInfo[uBi], "\t");
        uBOutTmp.getBlastInfo(uBInfoStr[0], uBInfoStr[1], str2int(uBInfoStr[3]), str2int(uBInfoStr[6]), str2int(uBInfoStr[7]), str2int(uBInfoStr[8]), str2int(uBInfoStr[9]), uBInfo[uBi]);
        uBOutTmp.contigOri = uBInfoStr[12];
        uBOut.push_back(uBOutTmp);
    }
    
    vector<string>().swap(uBInfo);
    vector<string>().swap(uBInfoStr);
    
    return uBOut;
}
