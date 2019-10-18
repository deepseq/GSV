/*
 
 FilterSummary.cpp
 
 This program is to do:
 1. summarize filtered results
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
*/

#include "FilterSummary.h"
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

////////////////////
/* Tool functions */
////////////////////

/* Confidence Estimation*/

float confEst(float cEConCov, float cEGenCov, float cEAlignLen, float cEAlignHit, float cEAlignDiff, float cEAlignGap)
{
    float bpCovRatio;

    if(cEGenCov == 0)
    {
        bpCovRatio = 10.0;
    }
    else
    {
        bpCovRatio = cEConCov/cEGenCov;
    }
    
    float alignLenLog = log10f(cEAlignLen);
    float alignHitExp = powf(10, cEAlignHit-2);
    float alignGapExp = powf(10, log10f((1-cEAlignGap/100)*cEAlignLen + 1));
    float alignDiff = cEAlignDiff + 1;
    
    float confValue = (bpCovRatio * alignLenLog)/(alignHitExp * alignGapExp * alignDiff) * 100;
    
    return confValue;
}


//////////////////////
/* filter functions */
//////////////////////

void filterSum(string fSID, vector<blastObj> &fSBlast, int fSDist)
{
    vector<blastObj>::iterator fSBlastFIter;
    vector<string> fSContigID;
    
    fSBlast = sortBlastContig(fSBlast);
    
    int fSi = 0;
    vector<blastObj>::iterator fSIter;
    vector<blastObj> fSBlastList;
    vector<blastObj> fSBlastListTmp;
    int fSListSize;
    
    regex fSRegStr;
    regex fSRegStrCov;
    
    vector<blastObj>::iterator fSAIter;
    vector<blastObj>::iterator fSBIter;
    vector<string> fSADataStr;
    vector<string> fSADataStrCov;
    vector<string> fSBDataStr;
    vector<string> fSBDataStrCov;
    
    float fSConfBPA;
    float fSConfBPB;
    
    // output file
    string fSFileStr;
    
    fSFileStr = fSID + "_GSVMiningRes/BlastRes/" + fSID + ".GRearr.candidates.summary.txt";
    ofstream fSFileBlastOut(fSFileStr);
    if(!fSFileBlastOut){printf("%s \n", "Error! Cannot open file to write! [filterSum]"); exit(1);}
    
    string fSFileOutStr;
    fSFileOutStr = "ID\tBP1_Chr\tBP1_GPos\tBP1_nonGPos\tBP2_Chr\tBP2_GPos\tBP2_nonGPos\tBP1_Cov\tBP1_Conf\tBP2_Cov\tBP2_Conf\tLength\tEvent\n";
    fSFileBlastOut << fSFileOutStr;
    fSFileOutStr.clear();
    
    // process
    for(fSBlastFIter = fSBlast.begin(); fSBlastFIter != fSBlast.end(); fSBlastFIter++)
    {
        fSContigID.push_back(fSBlastFIter->contigID);
    }
    
    fSContigID.erase(unique(fSContigID.begin(), fSContigID.end()), fSContigID.end());
    
    for(fSi = 0; fSi < int(fSContigID.size()); fSi++)
    {
        fSRegStr = fSContigID[fSi] + "\t";
        fSBlastList = grepBlastInfo(1, fSRegStr, fSBlast);
        fSListSize = int(fSBlastList.size());
        
        while(fSBlastList.size())
        {
            fSBlastListTmp.clear();
            fSBlastListTmp.push_back(*fSBlastList.begin());
            fSAIter = fSBlastListTmp.begin();
            
            fSBlastList.erase(fSBlastList.begin());
            
            for(fSIter = fSBlastList.begin(); fSIter != fSBlastList.end(); fSIter++)
            {
                fSBIter = fSIter;

                if(abs(fSAIter->contigPosB - fSBIter->contigPosA) < fSDist)
                {
                    fSRegStrCov = "--";
                    
                    fSADataStr = split(fSAIter->dataInfo, "\t");
                    if(regex_search(fSADataStr[4], fSRegStrCov))
                    {
                        fSADataStrCov = split(fSADataStr[4], "--");
                        fSADataStr[4] = fSADataStrCov[1];
                        
                        fSADataStrCov = split(fSADataStr[5], "--");
                        fSADataStr[5] = fSADataStrCov[1];
                    }
                    
                    fSBDataStr = split(fSBIter->dataInfo, "\t");
                    if(regex_search(fSBDataStr[4], fSRegStrCov))
                    {
                        fSBDataStrCov = split(fSBDataStr[4], "--");
                        fSBDataStr[4] = fSBDataStrCov[0];
                        
                        fSBDataStrCov = split(fSBDataStr[5], "--");
                        fSBDataStr[5] = fSBDataStrCov[0];
                    }
                    
                    fSConfBPA = confEst(str2float(fSADataStr[4]), str2float(fSADataStr[5]), str2float(fSADataStr[3]), float(fSListSize), float(abs(fSAIter->contigPosB - fSBIter->contigPosA)), str2float(fSADataStr[2]));
                    fSConfBPB = confEst(str2float(fSBDataStr[4]), str2float(fSBDataStr[5]), str2float(fSBDataStr[3]), float(fSListSize), float(abs(fSAIter->contigPosB - fSBIter->contigPosA)), str2float(fSBDataStr[2]));
                    
                    fSFileOutStr = fSAIter->contigID + "\t" + fSAIter->chrID + "\t" + fSADataStr[9] + "\t" + fSADataStr[8] + "\t" + fSBIter->chrID + "\t" + fSBDataStr[8] + "\t" + fSBDataStr[9] + "\t" + fSADataStr[4] + "\t" + float2str(fSConfBPA) + "\t" + fSBDataStr[4] + "\t" + float2str(fSConfBPB);

                    if(fSAIter->chrID == fSBIter->chrID)
                    {
                        fSFileOutStr = fSFileOutStr + "\t" + int2str(int(abs(fSAIter->chrPosB - fSBIter->chrPosA)));
                        
                        if(fSAIter->contigOri == fSBIter->contigOri)
                        {
                            if(fSAIter->contigOri == "forward")
                            {
                                if(fSAIter->chrPosB > fSBIter->chrPosA)
                                {
                                    fSFileOutStr = fSFileOutStr + "\t" + "Duplication" + "\n";
                                }
                                else
                                {
                                    fSFileOutStr = fSFileOutStr + "\t" + "Deletion" + "\n";
                                }
                            }
                            else
                            {
                                if(fSAIter->chrPosB > fSBIter->chrPosA)
                                {
                                    fSFileOutStr = fSFileOutStr + "\t" + "Deletion" + "\n";
                                }
                                else
                                {
                                    fSFileOutStr = fSFileOutStr + "\t" + "Duplication" + "\n";
                                }
                            }
                        }
                        else
                        {
                            fSFileOutStr = fSFileOutStr + "\t" + "Inversion" + "\n";
                        }
                    }
                    else
                    {
                        fSFileOutStr = fSFileOutStr + "\t" + "intrachr" + "\t" + "Translocation" + "\n";
                    }
                    
                    fSFileBlastOut << fSFileOutStr;
                    fSFileOutStr.clear();
                }
                else
                {
                    if(abs(fSAIter->chrPosB - fSBIter->chrPosA) < fSDist)
                    {
                        fSRegStrCov = "--";
                        
                        fSADataStr = split(fSAIter->dataInfo, "\t");
                        if(regex_search(fSADataStr[4], fSRegStrCov))
                        {
                            fSADataStrCov = split(fSADataStr[4], "--");
                            fSADataStr[4] = fSADataStrCov[1];
                            
                            fSADataStrCov = split(fSADataStr[5], "--");
                            fSADataStr[5] = fSADataStrCov[1];
                        }
                        
                        fSBDataStr = split(fSBIter->dataInfo, "\t");
                        if(regex_search(fSBDataStr[4], fSRegStrCov))
                        {
                            fSBDataStrCov = split(fSBDataStr[4], "--");
                            fSBDataStr[4] = fSBDataStrCov[0];
                            
                            fSBDataStrCov = split(fSBDataStr[5], "--");
                            fSBDataStr[5] = fSBDataStrCov[0];
                        }
                        
                        fSConfBPA = confEst(str2float(fSADataStr[4]), str2float(fSADataStr[5]), str2float(fSADataStr[3]), float(fSListSize), float(abs(fSAIter->contigPosB - fSBIter->contigPosA)), str2float(fSADataStr[2]));
                        fSConfBPB = confEst(str2float(fSBDataStr[4]), str2float(fSBDataStr[5]), str2float(fSBDataStr[3]), float(fSListSize), float(abs(fSAIter->contigPosB - fSBIter->contigPosA)), str2float(fSBDataStr[2]));
                        
                        fSFileOutStr = fSAIter->contigID + "\t" + fSAIter->chrID + "\t" + fSADataStr[9] + "\t" + fSADataStr[8] + "\t" + fSBIter->chrID + "\t" + fSBDataStr[8] + "\t" + fSBDataStr[9] + "\t" + fSADataStr[4] + "\t" + float2str(fSConfBPA) + "\t" + fSBDataStr[4] + "\t" + float2str(fSConfBPB);

                        fSFileOutStr = fSFileOutStr + "\t" + "ins" + "\t" + "Insertion" + "\n";
                        
                        fSFileBlastOut << fSFileOutStr;
                        fSFileOutStr.clear();
                    }
                }
            }
        }
        
        fSBlastList.clear();
        fSListSize = 0;
    }

    fSFileBlastOut.close();
}

////////////////////
/* main functions */
////////////////////

void GSVFilterSummary(string GSVFSID, int GSVFSDist)
{
    string GSVFSLine;
    string GSVFSFileStr;
    vector<string> GSVFSDataStr;
    
    blastObj GSVFSBlastObj;
    vector<blastObj> GSVFSBlastRes;
    GSVFSFileStr = GSVFSID + "_GSVMiningRes/BlastRes/" + GSVFSID + ".GRearr.candidates.bpCov.txt";
    ifstream GSVFSFileBlast(GSVFSFileStr);
    if(!GSVFSFileBlast){printf("%s \n", "Error! Cannot open blast result list file to read! [GSVFilterSum]"); exit(1);}
    while(1)
    {
        getline(GSVFSFileBlast, GSVFSLine); // careful! EOF check must be after this cmd as otherwise cannot get right EOF position.
        if(GSVFSFileBlast.eof()) break;
        GSVFSLine.erase(GSVFSLine.find_last_not_of(" \n\r\t")+1);
        
        GSVFSDataStr = split(GSVFSLine, "\t");
        GSVFSBlastObj.getBlastInfo(GSVFSDataStr[0], GSVFSDataStr[1], str2int(GSVFSDataStr[3]), str2int(GSVFSDataStr[6]), str2int(GSVFSDataStr[7]), str2int(GSVFSDataStr[8]), str2int(GSVFSDataStr[9]), GSVFSLine);
        GSVFSBlastObj.contigOri = GSVFSDataStr[12];
        GSVFSBlastRes.push_back(GSVFSBlastObj);
    }
    
    GSVFSFileBlast.close();
    
    filterSum(GSVFSID, GSVFSBlastRes, GSVFSDist);
}

