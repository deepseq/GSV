/*
 
 BlastFilter.cpp
 
 This program is to do:
 1. fitler false positive results
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
*/

#include "CurrentTime.h"
#include "CmdExecution.h"
#include "BlastFilter.h"
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

//////////////////////
/* filter functions */
//////////////////////

void rpkmFilter(vector<blastObj> &rFBlast, vector<assemblyObj> &rFAss, float rFThr)
{
    vector<blastObj> rFBlastF;
    vector<assemblyObj> rFAssF;
    
    int rFi;
    vector<blastObj> rFBlastTmp;
    regex rFRegStr;
    
    for(rFi = 0; rFi < int(rFAss.size()); rFi++)
    {
        if(rFAss[rFi].contigRPKM >= rFThr)
        {
            rFAssF.push_back(rFAss[rFi]);
            
            rFBlastTmp.clear();
            rFRegStr = rFAss[rFi].contigID + "\t";
            rFBlastTmp = grepBlastInfo(1, rFRegStr, rFBlast);
            rFBlastF.insert(rFBlastF.end(), rFBlastTmp.begin(), rFBlastTmp.end());
        }
    }
    
    rFBlast.clear();
    rFBlast = rFBlastF;
    
    rFAss.clear();
    rFAss = rFAssF;
    
    // free vector
    vector<blastObj>().swap(rFBlastF);
    vector<assemblyObj>().swap(rFAssF);
}

void spanFilter(vector<blastObj> &sFBlast, vector<assemblyObj> &sFAss, int sFThr)
{
    vector<blastObj> sFBlastFTmp;
    vector<blastObj> sFBlastF;
    vector<assemblyObj> sFAssF;
    
    int sFi;
    vector<blastObj>::iterator sFIter;
    vector<blastObj> sFBlastList;
    regex sFRegStr;
    vector<blastObj>::iterator sFAIter;
    vector<blastObj>::iterator sFBIter;
    vector<string> sFBlastInfoA;
    vector<string> sFBlastInfoB;
    vector<int> sFBlastConPos;
    vector<int> sFBlastChrPos;
    
    // blast results filtering
    for(sFi = 0; sFi < int(sFAss.size()); sFi++)
    {
        sFRegStr = sFAss[sFi].contigID + "\t";
        sFBlastList = grepBlastInfo(1, sFRegStr, sFBlast);
        
        // because the genomic position never repeats
        sFBlastList = sortBlastChr(sFBlastList);
        sFIter = sFBlastList.begin();
        while(sFIter != sFBlastList.end())
        {
            sFBIter = sFIter;
            sFBlastInfoB = split(sFIter->dataInfo, "\t");
            
            if(sFBlastInfoA.empty())
            {
                sFAIter = sFBIter;
                sFBlastInfoA = sFBlastInfoB;
                
                sFIter++;
                sFBIter = sFIter;
                sFBlastInfoB = split(sFIter->dataInfo, "\t");
            }
            
            sFBlastConPos.push_back(str2int(sFBlastInfoA[6]));
            sFBlastConPos.push_back(str2int(sFBlastInfoA[7]));
            sFBlastConPos.push_back(str2int(sFBlastInfoB[6]));
            sFBlastConPos.push_back(str2int(sFBlastInfoB[7]));
            sort(sFBlastConPos.begin(), sFBlastConPos.end());
            
            sFBlastChrPos.push_back(str2int(sFBlastInfoA[8]));
            sFBlastChrPos.push_back(str2int(sFBlastInfoA[9]));
            sFBlastChrPos.push_back(str2int(sFBlastInfoB[8]));
            sFBlastChrPos.push_back(str2int(sFBlastInfoB[9]));
            sort(sFBlastChrPos.begin(), sFBlastChrPos.end());
            
            if((sFBlastInfoA[1] == sFBlastInfoB[1])&&(sFBlastInfoA[12] == sFBlastInfoB[12]))
            {
                if((abs(sFBlastConPos[2] - sFBlastConPos[1]) < sFThr)&&(abs(sFBlastChrPos[2] - sFBlastChrPos[1]) < sFThr))
                {
                    sFBIter->contigAlign = str2int(sFBlastInfoA[3]) + str2int(sFBlastInfoB[3]);
                    sFBlastInfoB[3] = int2str(sFBIter->contigAlign);
                    
                    sFBIter->contigPosA = sFBlastConPos[0];
                    sFBlastInfoB[6] = int2str(sFBlastConPos[0]);
                    sFBIter->contigPosB = sFBlastConPos[3];
                    sFBlastInfoB[7] = int2str(sFBlastConPos[3]);
                    
                    if(sFBlastInfoA[12] == "forward")
                    {
                        sFBIter->chrPosA = sFBlastChrPos[0];
                        sFBlastInfoB[8] = int2str(sFBlastChrPos[0]);
                        sFBIter->chrPosB = sFBlastChrPos[3];
                        sFBlastInfoB[9] = int2str(sFBlastChrPos[3]);

                    }
                    else
                    {
                        sFBIter->chrPosA = sFBlastChrPos[3];
                        sFBlastInfoB[8] = int2str(sFBlastChrPos[3]);
                        sFBIter->chrPosB = sFBlastChrPos[0];
                        sFBlastInfoB[9] = int2str(sFBlastChrPos[0]);
                    }
                }
                else
                {
                    sFAIter->dataInfo = joinstring(sFBlastInfoA, "\t");
                    sFBlastFTmp.push_back(*sFAIter);
                }
            }
            else
            {
                sFAIter->dataInfo = joinstring(sFBlastInfoA, "\t");
                sFBlastFTmp.push_back(*sFAIter);
            }
            
            sFAIter = sFBIter;
            sFBlastInfoA = sFBlastInfoB;
            sFIter++;
            
            sFBlastConPos.clear();
            sFBlastChrPos.clear();
        }
        
        sFBIter->dataInfo = joinstring(sFBlastInfoB, "\t");
        sFBlastFTmp.push_back(*sFBIter);
        
        if(sFBlastFTmp.size() > 1)
        {
            sFBlastF.insert(sFBlastF.end(), sFBlastFTmp.begin(), sFBlastFTmp.end());
            sFAssF.push_back(sFAss[sFi]);
        }
        
        sFBlastFTmp.clear();
        sFBlastList.clear();
        sFBlastInfoA.clear();
        sFBlastInfoB.clear();
        sFBlastConPos.clear();
        sFBlastChrPos.clear();
    }
    
    sFBlast.clear();
    sFBlast = sFBlastF;
    
    sFAss.clear();
    sFAss = sFAssF;
    
    // free vector
    vector<blastObj>().swap(sFBlastList);
    vector<blastObj>().swap(sFBlastFTmp);
    vector<blastObj>().swap(sFBlastF);
    vector<assemblyObj>().swap(sFAssF);
}

void lengthFilter(vector<blastObj> &lFBlast, vector<assemblyObj> &lFAss, int lFThr)
{
    vector<blastObj> lFBlastF;
    vector<blastObj> lFBlastFTmp;
    vector<assemblyObj> lFAssF;

    int lFi = 0;
    vector<blastObj>::iterator lFIter;
    vector<blastObj> lFBlastList;
    regex lFRegStr;
    
    for(lFi = 0; lFi < int(lFAss.size()); lFi++)
    {
        lFRegStr = lFAss[lFi].contigID + "\t";
        lFBlastList = grepBlastInfo(1, lFRegStr, lFBlast);
        
        for(lFIter = lFBlastList.begin(); lFIter != lFBlastList.end(); lFIter++)
        {
            if(lFIter->contigAlign >= lFThr)
            {
                lFBlastFTmp.push_back(*lFIter);
            }
        }
        
        if(lFBlastFTmp.size() > 1)
        {
            lFBlastF.insert(lFBlastF.end(), lFBlastFTmp.begin(), lFBlastFTmp.end());
            lFAssF.push_back(lFAss[lFi]);
        }
        
        lFBlastFTmp.clear();
        lFBlastList.clear();
    }
    
    lFBlast.clear();
    lFBlast = lFBlastF;
    
    lFAss.clear();
    lFAss = lFAssF;
    
    // free vector
    vector<blastObj>().swap(lFBlastList);
    vector<blastObj>().swap(lFBlastFTmp);
    vector<blastObj>().swap(lFBlastF);
    vector<assemblyObj>().swap(lFAssF);
}

void gapFilter(vector<blastObj> &gFBlast, vector<assemblyObj> &gFAss, int gFThr)
{
    vector<blastObj> gFBlastFTmp;
    vector<blastObj> gFBlastF;
    vector<assemblyObj> gFAssF;
    
    int gFi;
    vector<blastObj>::iterator gFIter;
    vector<blastObj> gFBlastList;
    vector<blastObj> gFGapList;
    vector<blastObj> gFGapRes;
    regex gFRegStr;
    vector<blastObj>::iterator gFAIter;
    vector<blastObj>::iterator gFBIter;
    vector<string> gFBlastInfoA;
    vector<string> gFBlastInfoB;
    vector<int> gFBlastPos;
    
    int gFSign = 0;
    int gFGapSign = 0;
    
    for(gFi = 0; gFi < int(gFAss.size()); gFi++)
    {
        gFRegStr = gFAss[gFi].contigID + "\t";
        gFBlastList = grepBlastInfo(1, gFRegStr, gFBlast);
        
        gFBlastList = sortBlastContig(gFBlastList);
        gFIter = gFBlastList.begin();
        while(gFIter != gFBlastList.end())
        {
            gFBIter = gFIter;
            gFBlastInfoB = split(gFIter->dataInfo, "\t");
            
            if(gFBlastInfoA.empty())
            {
                gFAIter = gFBIter;
                gFBlastInfoA = gFBlastInfoB;
                
                gFIter++;
                gFBIter = gFIter;
                gFBlastInfoB = split(gFIter->dataInfo, "\t");
            }
            
            if((abs(gFAIter->contigPosA - gFBIter->contigPosA) < gFThr)||(abs(gFAIter->contigPosB - gFBIter->contigPosB) < gFThr))
            {
                gFGapList.push_back(*gFAIter);
                gFBlastPos.push_back(gFAIter->contigPosB);
                
                gFSign = 1;
            }
            else
            {
                gFGapList.push_back(*gFAIter);
                gFBlastPos.push_back(gFAIter->contigPosB);
                sort(gFBlastPos.begin(), gFBlastPos.end());
                
                if(abs(gFBlastPos.back() - gFBIter->contigPosA) < gFThr)
                {
                    gFGapRes.insert(gFGapRes.end(), gFGapList.begin(), gFGapList.end());
                    gFGapSign = 1;
                }
                else
                {
                    if(gFGapSign == 1)
                    {
                        gFGapRes.insert(gFGapRes.end(), gFGapList.begin(), gFGapList.end());
                        gFGapSign = 0;
                    }
                }
                
                gFGapList.clear();
                gFBlastPos.clear();
                gFSign = 0;
            }
            
            gFAIter = gFBIter;
            gFBlastInfoA = gFBlastInfoB;
            gFIter++;
        }
        
        if(gFSign == 1)
        {
            if(gFGapSign == 1)
            {
                gFGapRes.insert(gFGapRes.end(), gFGapList.begin(), gFGapList.end());
                gFGapRes.push_back(*gFBIter);
            }
        }
        else
        {
            if(gFGapSign == 1)
            {
                gFGapRes.push_back(*gFBIter);
            }
        }
         
        gFSign = 0;
        gFGapSign = 0;
        gFBlastPos.clear();
        gFGapList.clear();
        gFBlastInfoA.clear();
        gFBlastInfoB.clear();
        
        gFBlastList = sortBlastChr(gFBlastList);
        gFIter = gFBlastList.begin();
        while(gFIter != gFBlastList.end())
        {
            gFBIter = gFIter;
            gFBlastInfoB = split(gFIter->dataInfo, "\t");
            
            if(gFBlastInfoA.empty())
            {
                gFAIter = gFBIter;
                gFBlastInfoA = gFBlastInfoB;
                
                gFIter++;
                gFBIter = gFIter;
                gFBlastInfoB = split(gFIter->dataInfo, "\t");
            }
            
            if((abs(gFAIter->chrPosA - gFBIter->chrPosA) < gFThr)||(abs(gFAIter->chrPosB - gFBIter->chrPosB) < gFThr))
            {
                gFGapList.push_back(*gFAIter);
                gFBlastPos.push_back(gFAIter->chrPosB);
                
                gFSign = 1;
            }
            else
            {
                gFGapList.push_back(*gFAIter);
                gFBlastPos.push_back(gFAIter->chrPosB);
                sort(gFBlastPos.begin(), gFBlastPos.end());
                
                if(abs(gFBlastPos.back() - gFBIter->chrPosA) < gFThr)
                {
                    gFGapRes.insert(gFGapRes.end(), gFGapList.begin(), gFGapList.end());
                    gFGapSign = 1;
                }
                else
                {
                    if(gFGapSign == 1)
                    {
                        gFGapRes.insert(gFGapRes.end(), gFGapList.begin(), gFGapList.end());
                        gFGapSign = 0;
                    }
                }
                
                gFGapList.clear();
                gFBlastPos.clear();
                gFSign = 0;
            }
            
            gFAIter = gFBIter;
            gFBlastInfoA = gFBlastInfoB;
            gFIter++;
        }
        
        if(gFSign == 1)
        {
            if(gFGapSign == 1)
            {
                gFGapRes.insert(gFGapRes.end(), gFGapList.begin(), gFGapList.end());
                gFGapRes.push_back(*gFBIter);
            }
        }
        else
        {
            if(gFGapSign == 1)
            {
                gFGapRes.push_back(*gFBIter);
            }
        }
        
        gFSign = 0;
        gFGapSign = 0;
        gFBlastPos.clear();
        gFGapList.clear();
        gFBlastInfoA.clear();
        gFBlastInfoB.clear();

        gFGapRes = sortBlastContig(gFGapRes);
        gFBlastFTmp = uniqueBlast(gFGapRes);
        
        gFBlastF.insert(gFBlastF.end(), gFBlastFTmp.begin(), gFBlastFTmp.end());
         
        gFGapRes.clear();
        gFBlastList.clear();
        gFBlastFTmp.clear();
    }
    
    for(gFi = 0; gFi < int(gFAss.size()); gFi++)
    {
        gFRegStr = gFAss[gFi].contigID + "\t";
        gFBlastList = grepBlastInfo(1, gFRegStr, gFBlastF);
        
        if(!gFBlastList.empty())
        {
            gFAssF.push_back(gFAss[gFi]);
        }
    }
    
    gFBlast.clear();
    gFBlast = gFBlastF;
    
    gFAss.clear();
    gFAss = gFAssF;
    
    vector<blastObj>().swap(gFGapRes);
    vector<blastObj>().swap(gFGapList);
    vector<blastObj>().swap(gFBlastList);
    vector<blastObj>().swap(gFBlastFTmp);
    vector<blastObj>().swap(gFBlastF);
    vector<assemblyObj>().swap(gFAssF);
}

////////////////////
/* main functions */
////////////////////

void GSVBlastFilter(string GSVBFID, int GSVBFGap, int GSVBFDist)
{
    // read files
    
    string GSVBFLine;
    string GSVBFFileStr;
    vector<string> GSVBFDataStr;
    regex GSVBFRegStr;
    
    GSVBFRegStr = "contig_";
    
    vector<string> GSVBFContigID;
    blastObj GSVBFBlastObj;
    vector<blastObj> GSVBFBlastRes;
    GSVBFFileStr = GSVBFID + "_GSVMiningRes/BlastRes/" + GSVBFID + ".contig.blast.txt";
    ifstream GSVBFFileBlast(GSVBFFileStr);
    if(!GSVBFFileBlast){printf("%s \n", "Error! Cannot open file to read! [GSVBlastFilter]"); exit(1);}
    while(1)
    {
        getline(GSVBFFileBlast, GSVBFLine); 
        if(GSVBFFileBlast.eof()) break;
        GSVBFLine.erase(GSVBFLine.find_last_not_of(" \n\r\t")+1);
        
        if(regex_search(GSVBFLine, GSVBFRegStr, regex_constants::match_continuous))
        {
            GSVBFDataStr = split(GSVBFLine, "\t");
            GSVBFBlastObj.getBlastInfo(GSVBFDataStr[0], GSVBFDataStr[1], str2int(GSVBFDataStr[3]), str2int(GSVBFDataStr[6]), str2int(GSVBFDataStr[7]), str2int(GSVBFDataStr[8]), str2int(GSVBFDataStr[9]), GSVBFLine);
            GSVBFBlastRes.push_back(GSVBFBlastObj);
            GSVBFContigID.push_back(GSVBFDataStr[0]);
        }
    }
    
    GSVBFFileBlast.close();
    
    GSVBFRegStr = ">contig_";
    
    assemblyObj GSVBFAssObj;
    vector<assemblyObj> GSVBFAssRes;
    GSVBFFileStr = GSVBFID + "_GSVMiningRes/AssemblyRes/" + GSVBFID + ".contig.fasta";
    ifstream GSVBFFileAss(GSVBFFileStr);
    if(!GSVBFFileAss){printf("%s \n", "Error! Cannot open file to read! [GSVBlastFilter]"); exit(1);}
    
    while(1)
    {
        getline(GSVBFFileAss, GSVBFLine); // careful! EOF check must be after this cmd as otherwise cannot get right EOF position.
        if(GSVBFFileAss.eof()) break;
        GSVBFLine.erase(GSVBFLine.find_last_not_of(" \n\r\t")+1);
        
        if(regex_search(GSVBFLine, GSVBFRegStr, regex_constants::match_continuous))
        {
            GSVBFLine.erase(0,1);
            GSVBFAssObj.contigID = GSVBFLine;
        }
        else
        {
            GSVBFAssObj.contigSeq = GSVBFLine;
            GSVBFAssObj.contigLen = int(GSVBFLine.length());
            GSVBFAssRes.push_back(GSVBFAssObj);
        }
    }
    
    GSVBFFileAss.close();
    
    vector<string>().swap(GSVBFDataStr);

    GSVBFFileStr = GSVBFID + "_GSVMiningRes/BlastRes/" + GSVBFID + ".GRearr.candidates.txt";
    ofstream GSVBFBlastOut(GSVBFFileStr);
    if(!GSVBFBlastOut){printf("%s \n", "Error! Cannot open file to write! [GSVBlastFilter]"); exit(1);}
    
    GSVBFFileStr = GSVBFID + "_GSVMiningRes/AssemblyRes/" + GSVBFID + ".GRearr.candidates.fasta";
    ofstream GSVBFFastaOut(GSVBFFileStr);
    if(!GSVBFFastaOut){printf("%s \n", "Error! Cannot open file to write! [GSVBlastFilter]"); exit(1);}
    
    GSVBFContigID.erase(unique(GSVBFContigID.begin(), GSVBFContigID.end()), GSVBFContigID.end());
    
    int GSVBFi;
    int GSVBFj;
    int GSVBFBlastLen = 0;
    int GSVBFContigLen = 0;
    vector<string>::iterator GSVBFIter;
    vector<blastObj> GSVBFBlastTmp;
    for(GSVBFIter = GSVBFContigID.begin(); GSVBFIter != GSVBFContigID.end(); GSVBFIter++)
    {
        GSVBFBlastTmp.clear();

        GSVBFRegStr = *GSVBFIter + "\t";
        GSVBFBlastTmp = grepBlastInfo(1, GSVBFRegStr, GSVBFBlastRes);
        
        if(GSVBFBlastTmp.size() > 1)
        {
            for(GSVBFi = 0; GSVBFi < int(GSVBFAssRes.size()); GSVBFi++)
            {
                if(GSVBFAssRes[GSVBFi].contigID == *GSVBFIter)
                {
                    GSVBFContigLen = GSVBFAssRes[GSVBFi].contigLen;
                    GSVBFFastaOut << ">" << GSVBFAssRes[GSVBFi].contigID << "\n";
                    GSVBFFastaOut << GSVBFAssRes[GSVBFi].contigSeq << "\n";
                    break;
                }
            }
            
            GSVBFBlastTmp = sortBlastContig(GSVBFBlastTmp);
            
            for(GSVBFj = 0; GSVBFj < (int(GSVBFBlastTmp.size())-1); GSVBFj++)
            {
                if(abs(GSVBFBlastTmp[GSVBFj+1].contigPosA - GSVBFBlastTmp[GSVBFj].contigPosB) > GSVBFGap)
                {
                    break;
                }
            }
            
            for(GSVBFi = 0; GSVBFi < int(GSVBFBlastTmp.size()); GSVBFi++)
            {
                GSVBFBlastLen = GSVBFBlastLen + GSVBFBlastTmp[GSVBFi].contigAlign;
            }
            
            for(GSVBFi = 0; GSVBFi < int(GSVBFBlastTmp.size()); GSVBFi++)
            {
                GSVBFBlastOut << GSVBFBlastTmp[GSVBFi].dataInfo;
                
                if((GSVBFBlastTmp[GSVBFi].chrPosA - GSVBFBlastTmp[GSVBFi].chrPosB) > 0)
                {
                    GSVBFBlastOut << "\t" << "reverse";
                }
                else
                {
                    GSVBFBlastOut << "\t" << "forward";
                }
                
                if(GSVBFj == (int(GSVBFBlastTmp.size())-1))
                {
                    if(abs(GSVBFBlastLen - GSVBFContigLen) > GSVBFDist)
                    {
                        GSVBFBlastOut << "\t" << "partially_aligned";
                    }
                    else
                    {
                        GSVBFBlastOut << "\t" << "wholly_aligned";
                    }
                }
                else
                {
                    GSVBFBlastOut << "\t" << "partially_aligned";
                }
                
                GSVBFBlastOut << "\n";
            }
        }
        
        GSVBFBlastLen = 0;
    }
    
    vector<blastObj>().swap(GSVBFBlastTmp);
    
    GSVBFFastaOut.close();
    GSVBFBlastOut.close();
}

void GSVPrimaryFilter(string GSVPFID, int GSVPFGap, int GSVPFDist, float GSVPFRPKM, int GSVPFSpan)
{
    // read files
    
    string GSVPFTmp;
    char *GSVPFCmd = new char[1024];
    
    string GSVPFLine;
    string GSVPFFileStr;
    vector<string> GSVPFDataStr;
    
    blastObj GSVPFBlastObj;
    vector<blastObj> GSVPFBlastRes;
    GSVPFFileStr = GSVPFID + "_GSVMiningRes/BlastRes/" + GSVPFID + ".GRearr.candidates.txt";
    ifstream GSVPFFileBlast(GSVPFFileStr);
    if(!GSVPFFileBlast){printf("%s \n", "Error! Cannot open blast result list file to read! [GSVPrimaryFilter]"); exit(1);}
    while(1)
    {
        getline(GSVPFFileBlast, GSVPFLine); // careful! EOF check must be after this cmd as otherwise cannot get right EOF position.
        if(GSVPFFileBlast.eof()) break;
        GSVPFLine.erase(GSVPFLine.find_last_not_of(" \n\r\t")+1);
        
        GSVPFDataStr = split(GSVPFLine, "\t");
        GSVPFBlastObj.getBlastInfo(GSVPFDataStr[0], GSVPFDataStr[1], str2int(GSVPFDataStr[3]), str2int(GSVPFDataStr[6]), str2int(GSVPFDataStr[7]), str2int(GSVPFDataStr[8]), str2int(GSVPFDataStr[9]), GSVPFLine);
        GSVPFBlastObj.contigOri = GSVPFDataStr[12];
        GSVPFBlastRes.push_back(GSVPFBlastObj);
    }
    
    GSVPFFileBlast.close();
    
    assemblyObj GSVPFAssObj;
    vector<assemblyObj> GSVPFAssRes;
    GSVPFFileStr = GSVPFID + "_GSVMiningRes/ReMappingRes/" + GSVPFID + ".contigs.rpkm.txt";
    ifstream GSVPFFileAss(GSVPFFileStr);
    if(!GSVPFFileAss){printf("%s \n", "Error! Cannot open contig remapping rpkm file to read! [GSVPrimaryFilter]"); exit(1);}
    
    while(1)
    {
        getline(GSVPFFileAss, GSVPFLine); // careful! EOF check must be after this cmd as otherwise cannot get right EOF position.
        if(GSVPFFileAss.eof()) break;
        GSVPFLine.erase(GSVPFLine.find_last_not_of(" \n\r\t")+1);
        
        GSVPFDataStr = split(GSVPFLine, "\t");
        GSVPFAssObj.getAssInfo(GSVPFDataStr[0], str2int(GSVPFDataStr[1]), str2longint(GSVPFDataStr[2]), str2float(GSVPFDataStr[3]), str2float(GSVPFDataStr[4]));
        GSVPFAssRes.push_back(GSVPFAssObj);
    }
    
    GSVPFFileAss.close();
    
    vector<string>().swap(GSVPFDataStr);

    cout << ">> " << currentDateTime().c_str() << " <<\t--------> RPKM Filtering\n";
    rpkmFilter(GSVPFBlastRes, GSVPFAssRes, GSVPFRPKM);
    
    cout << ">> " << currentDateTime().c_str() << " <<\t--------> SPAN Filtering\n";
    spanFilter(GSVPFBlastRes, GSVPFAssRes, GSVPFSpan);
    
    cout << ">> " << currentDateTime().c_str() << " <<\t--------> LEN Filtering\n";
    lengthFilter(GSVPFBlastRes, GSVPFAssRes, GSVPFGap);
    
    cout << ">> " << currentDateTime().c_str() << " <<\t--------> GAP Filtering\n";
    gapFilter(GSVPFBlastRes, GSVPFAssRes, GSVPFDist);
    
    int GSVPFi;
    
    string GSVPFBOut = GSVPFID + "_GSVMiningRes/BlastRes/" + GSVPFID + ".GRearr.candidates.filter.txt";
    ofstream GSVPFBlastOut(GSVPFBOut);
    if(!GSVPFBlastOut){printf("%s \n", "Error! Cannot open blast result list file to write! [GSVPrimaryFilter]"); exit(1);}
    
    for(GSVPFi = 0; GSVPFi < int(GSVPFBlastRes.size()); GSVPFi++)
    {
        GSVPFBlastOut << GSVPFBlastRes[GSVPFi].dataInfo << "\n";
    }
    
    GSVPFBlastOut.close();
    
    string GSVPFAOut = GSVPFID + "_GSVMiningRes/BlastRes/" + GSVPFID + ".GRearr.candidates.filter.rpkm.txt";
    ofstream GSVPFAssOut(GSVPFAOut);
    if(!GSVPFAssOut){printf("%s \n", "Error! Cannot open blast result rpkm file to write! [GSVPrimaryFilter]"); exit(1);}
    
    for(GSVPFi = 0; GSVPFi < int(GSVPFAssRes.size()); GSVPFi++)
    {
        GSVPFAssOut << GSVPFAssRes[GSVPFi].contigID << "\t";
        GSVPFAssOut << GSVPFAssRes[GSVPFi].contigLen << "\t";
        GSVPFAssOut << GSVPFAssRes[GSVPFi].contigCount << "\t";
        GSVPFAssOut << GSVPFAssRes[GSVPFi].contigRPKM << "\t";
        GSVPFAssOut << GSVPFAssRes[GSVPFi].contigLogRPKM << "\n";
    }
    
    GSVPFAssOut.close();
    
    cout << ">> " << currentDateTime().c_str() << " <<\t--------> ALIGN Filtering\n";
    
    GSVPFTmp = "awk '{print $1}' " + GSVPFID + "_GSVMiningRes/BlastRes/" + GSVPFID + ".GRearr.candidates.filter.txt | uniq -c | awk '{if($1 <= 8){print $2}}' > " + GSVPFID + "_GSVMiningRes/BlastRes/" + GSVPFID + ".GRearr.candidates.filter.id.txt";
    strcpy(GSVPFCmd, GSVPFTmp.c_str());
    cmdExec(GSVPFCmd, 0);
    
    GSVPFTmp = "awk 'NR==FNR{a[$1]=$1;next}{if(a[$1]){print $0}}' " + GSVPFID + "_GSVMiningRes/BlastRes/" + GSVPFID + ".GRearr.candidates.filter.id.txt " + GSVPFID + "_GSVMiningRes/BlastRes/" + GSVPFID + ".GRearr.candidates.filter.txt > " + GSVPFID + "_GSVMiningRes/BlastRes/" + GSVPFID + ".GRearr.candidates.filter.rm.txt";
    strcpy(GSVPFCmd, GSVPFTmp.c_str());
    cmdExec(GSVPFCmd, 0);
    
    GSVPFTmp = "awk 'NR==FNR{a[$1]=$1;next}{if(a[$1]){print $0}}' " + GSVPFID + "_GSVMiningRes/BlastRes/" + GSVPFID + ".GRearr.candidates.filter.id.txt " + GSVPFID + "_GSVMiningRes/BlastRes/" + GSVPFID + ".GRearr.candidates.filter.rpkm.txt > " + GSVPFID + "_GSVMiningRes/BlastRes/" + GSVPFID + ".GRearr.candidates.filter.rm.rpkm.txt";
    strcpy(GSVPFCmd, GSVPFTmp.c_str());
    cmdExec(GSVPFCmd, 0);
    
    GSVPFTmp = "rm " + GSVPFID + "_GSVMiningRes/BlastRes/" + GSVPFID + ".GRearr.candidates.filter.id.txt";
    strcpy(GSVPFCmd, GSVPFTmp.c_str());
    cmdExec(GSVPFCmd, 0);
    
    delete[] GSVPFCmd;
}
