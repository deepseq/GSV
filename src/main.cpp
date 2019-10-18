/*
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
*/

////////////////
/* head files */
////////////////

#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>

#include "CurrentTime.h"
#include "CmdExecution.h"
#include "StringOperation.h"
#include "DataPreparation.h"
#include "MappingReads.h"
#include "BaseCoverage.h"
#include "MappingSummary.h"
#include "ReadsAssembly.h"
#include "BlastFilter.h"
#include "CoverageFilter.h"
#include "FilterSummary.h"

//////////////////////
/* global variables */
//////////////////////

const static char* helpInfo = "usage:\n"
                              "\tGSVMining -option variable ...\n"
                              "\nOptions:\n"
                              "\tRequired:\n"
                              "\t-i/--analysisIDList\tanalysis id list in text file (indentical to id of raw read file)\n"
                              "\t-c/--chromosomeList\tchromosome list in txt file (per chromosome per line)\n"
                              "\t-r/--readsDirectory\tdirectory of raw read files (ex: id.read_1.fastq id.read_2.fastq)\n"
                              "\t-p/--avgReadLength\taverage read length\n"
                              "\t-q/--avgInsertSize\taverage insert size\n"
                              "\t-b/--sdInsertSize\tstandard deviation of insert size\n"
                              "\t-f/--referencePath\tpath of reference in fasta format\n"
                              "\t-o/--outputDirectory\tdirectory of output\n"
                              "\tOptional:\n"
                              "\t-g/--gapThreshold\tgap allowed between two blast alignment [default: 12]\n"
                              "\t-d/--distThreshold\tdistance allowed between two blast alignment [default: 24]\n"
                              "\t-e/--rpkmThreshold\trpkm threshold to filter low coverage candidates [default: 1]\n"
                              "\t-s/--spanThreshold\tspan threshold to filter false positive candidates [default: 200]\n"
                              "\t-l/--maxFragLength\tmaximal fragment length [default: 1000]\n"
                              "\t-t/--numberThreads\tthe number of threads used in bwa, samtools [default: 12]\n"
                              "\t-k/--keepTmp\tkeep the temp files [default: no]";

const static char* noteInfo = "Note:\tPipeline for Genome Structural Variant Analysis\n"
                              "\t1. samtools, bedtools, bwa, assembly-tool, blast, bam2fastx, perl should be already installed.\n"
                              "\t2. blastdb should be built first.\n"
                              "\t3. reference fasta sequence should be indexed first (directory should be the same with blastdb).\n"
                              "\t4. read file should be named as ID.read_1.fastq.gz & ID.read_2.fastq.gz.\n";

const static char* authorInfo = "\tFei Sang\n"
                                "\tDeep Seq, Queen's Medical School\n"
                                "\tSchool of Life Sciences, University of Nottingham\n"
                                "\tNottingham, Nottinghamshire, NG7 2JE, United Kingdom\n"
                                "\tEmail: fei.sang@nottingham.ac.uk\n";

using namespace std;

const string versionInfo = "V1.4";

///////////////////
/* program class */
///////////////////

class GSVObj
{
public:
    string runID;
    string genomeRef;
    vector<string> chrList;
    
public:
    void getInfo(string subID, string subRef, vector<string> subList);
};

void GSVObj::getInfo(string subID, string subRef, vector<string> subList)
{
    runID = subID;
    genomeRef = subRef;
    chrList = subList;
}

////////////////////
/* main functions */
////////////////////

int main(int argNum, char *argVar[])
{
    ///////////////////
    /* Welcome Words */
    ///////////////////
    
    string welcomeWords = "Welcome to use GSVMining " + versionInfo;
    cout << string(welcomeWords.length()+20, '*') << "\n";
    cout << "*" << string(9, ' ') << welcomeWords << string(9, ' ') << "*" << "\n";
    cout << "*" << string(welcomeWords.length()+18, ' ') << "*" << "\n";
    cout << "*" << string(13, ' ') << "Fei Sang,DeepSeq,QMC" << string(14, ' ') << "*" << "\n";
    cout << "*" << string(7, ' ') << "Life Sciences,Univ. of Nottingham" << string(7, ' ') << "*" << "\n";
    cout << "*" << string(12, ' ') << "Nottingham, NG7 2UH, UK" << string(12, ' ') << "*" << "\n";
    cout << "*" << string(11, ' ') << "fei.sang@nottingham.ac.uk" << string(11, ' ') << "*" << "\n";
    cout << string(welcomeWords.length()+20, '*') << "\n";
    cout << "\n";
    
    ////////////////////////////////////////////
    /* 1. obtain variables from screen inputs */
    ////////////////////////////////////////////
    
    // claim variables
    int getOptIndex = 0;
    int getOptInfo;

    string runID;
    string chrList;
    string readsDir;
    string avgLength;
    string avgInsert;
    string sdInsert;
    string refPath;
    string outDir;
    
    int GSVGap = 12;  // also be used as length threshold to filter false positive
    int GSVDist = 24; // also be used as gap threshold to filter false positive
    float GSVRPKM = 1.00;
    int GSVSpan = 200;
    int GSVMaxLen = 1000;
    int GSVThreads = 12;
    string GSVKeepTmp = "no";

    static struct option longOptions[] =
    {
        {"help",            no_argument,       0, 'h'},
        {"version",         no_argument,       0, 'v'},
        {"analysisIDList",  required_argument, 0, 'i'},
        {"chromosomeList",  required_argument, 0, 'c'},
        {"readsDirectory",  required_argument, 0, 'r'},
        {"avgReadLength",   required_argument, 0, 'p'},
        {"avgInsertSize",   required_argument, 0, 'q'},
        {"sdInsertSize",    required_argument, 0, 'b'},
        {"referencePath",   required_argument, 0, 'f'},
        {"outputDirectory", required_argument, 0, 'o'},
        {"gapThreshold",    optional_argument, 0, 'g'},
        {"distThreshold",   optional_argument, 0, 'd'},
        {"rpkmThreshold",   optional_argument, 0, 'e'},
        {"spanThreshold",   optional_argument, 0, 's'},
        {"maxFragLength",   optional_argument, 0, 'l'},
        {"numberThreads",   optional_argument, 0, 't'},
        {"keepTmp",         optional_argument, 0, 'k'}
    };

    while(1)
    {
        // getopt_long stores the option index here
        int optionIndex = 0;
        getOptInfo = getopt_long(argNum, argVar, "hvi:c:r:p:q:b:f:o:g:d:e:s:l:t:k:", longOptions, &optionIndex);

        // detect the end of the options
        if(getOptInfo == -1)
        {
            printf("%s \n", helpInfo);
            printf("%s \n", noteInfo);
            printf("%s \n", authorInfo);
            break;
        }

        switch(getOptInfo)
        {
            case 'h':
            {
                printf("%s \n", helpInfo);
                printf("%s \n", noteInfo);
                printf("%s \n", authorInfo);
                exit(0);
            }

            case 'v':
            {
                printf("%s \n", "GSVMining Version: 0.1\n");
                exit(0);
            }

            case 'i':{runID     = optarg; getOptIndex++; break;}
            case 'c':{chrList   = optarg; getOptIndex++; break;}
            case 'r':{readsDir  = optarg; getOptIndex++; break;}
            case 'p':{avgLength = optarg; getOptIndex++; break;}
            case 'q':{avgInsert = optarg; getOptIndex++; break;}
            case 'b':{sdInsert  = optarg; getOptIndex++; break;}
            case 'f':{refPath   = optarg; getOptIndex++; break;}
            case 'o':{outDir    = optarg; getOptIndex++; break;}
                
            case 'g':
            {
                if(optarg != NULL){GSVGap = atoi(optarg);}
                getOptIndex++;
                break;
            }
                
            case 'd':
            {
                if(optarg != NULL){GSVDist = atoi(optarg);}
                getOptIndex++;
                break;
            }
                
            case 'e':
            {
                if(optarg != NULL){GSVRPKM = atof(optarg);}
                getOptIndex++;
                break;
            }
                
            case 's':
            {
                if(optarg != NULL){GSVSpan = atoi(optarg);}
                getOptIndex++;
                break;
            }
                
            case 'l':
            {
                if(optarg != NULL){GSVMaxLen = atoi(optarg);}
                getOptIndex++;
                break;
            }
                
            case 't':
            {
                if(optarg != NULL){GSVThreads = atoi(optarg);}
                getOptIndex++;
                break;
            }
                
            case 'k':
            {
                if(optarg != NULL){GSVKeepTmp = atoi(optarg);}
                getOptIndex++;
                break;
            }

            case '?':
            {
                // getopt_long already printed an error message.
                printf("%s \n\n", "Error!!! Please check parameters you provided!");
                return EXIT_FAILURE;
            }

            default: {abort();}
        }
    }

    // check if each variable has a value
    if(runID.empty())     {printf("%s \n\n", "Error!!! Please provide analysis ID list! [-i] or try [-h] to get help");                     return EXIT_FAILURE;}
    if(chrList.empty())   {printf("%s \n\n", "Error!!! Please provide chromosome list! [-c] or try [-h] to get help");                      return EXIT_FAILURE;}
    if(readsDir.empty())  {printf("%s \n\n", "Error!!! Please provide reads directory! [-r] or try [-h] to get help");                      return EXIT_FAILURE;}
    if(avgLength.empty()) {printf("%s \n\n", "Error!!! Please provide read length! [-p] or try [-h] to get help");                          return EXIT_FAILURE;}
    if(avgInsert.empty()) {printf("%s \n\n", "Error!!! Please provide insert size! [-q] or try [-h] to get help");                          return EXIT_FAILURE;}
    if(sdInsert.empty())  {printf("%s \n\n", "Error!!! Please provide standard deviation of insert size! [-b] or try [-h] to get help");    return EXIT_FAILURE;}
    if(refPath.empty())   {printf("%s \n\n", "Error!!! Please provide reference path! [-f] or try [-h] to get help");                       return EXIT_FAILURE;}
    if(outDir.empty())    {printf("%s \n\n", "Error!!! Please provide analysis output directory! [-o] or try [-h] to get help");            return EXIT_FAILURE;}

    // check index file exists
    string fileCheck = "";
    fileCheck = refPath + ".sa";
    
    struct stat buf;
    if(stat(fileCheck.c_str(), &buf) != 0)
    {
        printf("%s \n\n", "Error!!! fasta file is not indexed by bwa!");
        return EXIT_FAILURE;
    }
    
    fileCheck = refPath + ".nhr";
    if(stat(fileCheck.c_str(), &buf) != 0)
    {
        printf("%s \n\n", "Error!!! fasta file is not indexed by makeblastdb!");
        return EXIT_FAILURE;
    }
    
    ///////////////////
    /* 2. read files */
    ///////////////////
    
    cout << ">> " << currentDateTime().c_str() << " <<\tReading Files\n";
    
    // claim variables
    string sysTmp = "";
    vector<string> arrRunID;
    vector<string> arrChrList;

    ifstream fileRunID(runID.c_str());
    if(!fileRunID){printf("%s \n", "Error! Cannot open analysis ID list file to read!"); exit(1);}
    
    sysTmp = "";
    while(1)
    {
        getline(fileRunID, sysTmp); // careful! EOF check must be after this cmd as otherwise cannot get right EOF position.
        if(fileRunID.eof()) break;

        sysTmp.erase(sysTmp.find_last_not_of(" \n\r\t")+1);
        if(sysTmp != "") arrRunID.push_back(sysTmp);
    }
    fileRunID.close();
    
    ifstream fileChrList(chrList.c_str());
    if(!fileChrList){printf("%s \n", "Error! Cannot open chromosome list file to read!"); exit(1);}
    
    sysTmp = "";
    while(1)
    {
        getline(fileChrList, sysTmp); // careful! EOF check must be after this cmd as otherwise cannot get right EOF position.
        if(fileChrList.eof()) break;

        sysTmp.erase(sysTmp.find_last_not_of(" \n\r\t")+1);
        if(sysTmp != "") arrChrList.push_back(sysTmp);
    }
    fileChrList.close();

    /////////////////////
    /* 3. prepare data */
    /////////////////////
    
    cout << ">> " << currentDateTime().c_str() << " <<\tPreparing Data\n";
    
    // claim variables
    vector<string>::iterator runIter;
    
    string mainTmp;
    char *mainCmd = new char[1024];

    /////////////////////
    /* 4. analyze data */
    /////////////////////
    
    // claim variables
    
    long int totalMapped = 0;
    long int totalMappedContig = 0;
    long int totalMappedChr = 0;
    
    cout << ">> " << currentDateTime().c_str() << " <<\tAnalyzing Data\n\n";

    for(runIter = arrRunID.begin(); runIter != arrRunID.end(); runIter++)
    {
        GSVObj GSVTmp;
        GSVTmp.getInfo(*runIter, refPath, arrChrList);
        
        // create log files
        
        string GSVLogStr = GSVTmp.runID + ".log.txt";
        ofstream GSVLogOut(GSVLogStr);
        if(!GSVLogOut){printf("%s \n", "Error! Cannot open log file to write! Please check authority control."); exit(1);}
        
        // pipeline begins
        
        GSVLogOut << ">> " << currentDateTime().c_str() << " <<\tAnalyzing ID: " << GSVTmp.runID << "\n";
        cout << ">> " << currentDateTime().c_str() << " <<\tAnalyzing ID: " << GSVTmp.runID << "\n";
        dataPre(GSVTmp.runID, GSVTmp.genomeRef);

        GSVLogOut << ">> " << currentDateTime().c_str() << " <<\t----> 1. BWA Mapping & Filtering\n";
        cout << ">> " << currentDateTime().c_str() << " <<\t----> 1. BWA Mapping & Filtering\n";
        mapReads(refPath, readsDir, GSVTmp.runID, GSVThreads);
        
        GSVLogOut << ">> " << currentDateTime().c_str() << " <<\t----> 2. Counting Base Coverage for Each Chromosome\n";
        cout << ">> " << currentDateTime().c_str() << " <<\t----> 2. Counting Base Coverage for Each Chromosome\n";
        baseCov(GSVTmp.runID, GSVTmp.chrList);
        
        GSVLogOut << ">> " << currentDateTime().c_str() << " <<\t----> 3. Summerizing Mapping Results\n";
        cout << ">> " << currentDateTime().c_str() << " <<\t----> 3. Summerizing Mapping Results\n";
        totalMappedChr = mapSum(GSVTmp.runID);
        
        GSVLogOut << ">> " << currentDateTime().c_str() << " <<\t----> 4. De Novo Assembly\n";
        cout << ">> " << currentDateTime().c_str() << " <<\t----> 4. De Novo Assembly\n";
        totalMapped = readsAss(avgLength, avgInsert, sdInsert, refPath, GSVTmp.runID, totalMappedChr, GSVGap, GSVDist, GSVThreads);
        totalMappedContig = totalMapped - totalMappedChr;
        
        GSVLogOut << ">> " << currentDateTime().c_str() << " <<\t----> 5. Candidate Primary Filtering\n";
        cout << ">> " << currentDateTime().c_str() << " <<\t----> 5. Candidate Primary Filtering\n";
        GSVPrimaryFilter(GSVTmp.runID, GSVGap, GSVDist, GSVRPKM, GSVSpan);
        
        GSVLogOut << ">> " << currentDateTime().c_str() << " <<\t----> 6. Candidate Secondary Filtering\n";
        cout << ">> " << currentDateTime().c_str() << " <<\t----> 6. Candidate Secondary Filtering\n";
        GSVCovFilter(GSVTmp.runID, GSVDist, GSVMaxLen, totalMappedChr, totalMappedContig);

        GSVLogOut << ">> " << currentDateTime().c_str() << " <<\t----> 7. Candidate Summary\n";
        cout << ">> " << currentDateTime().c_str() << " <<\t----> 7. Candidate Summary\n";
        GSVFilterSummary(GSVTmp.runID, GSVDist);

        if(GSVKeepTmp == "yes")
        {
            mainTmp = "mv " + GSVTmp.runID + "_GSVMiningRes" + " " + outDir;
            strcpy(mainCmd, mainTmp.c_str());
            cmdExec(mainCmd, 0);
        }
        else
        {
            mainTmp = "mv " + GSVTmp.runID + "_GSVMiningRes/BlastRes/" + GSVTmp.runID + ".GRearr.candidates.summary.txt" + " " + outDir;
            strcpy(mainCmd, mainTmp.c_str());
            cmdExec(mainCmd, 0);
            
            mainTmp = "mv " + GSVTmp.runID + "_GSVMiningRes/MappingRes/*.basecov.txt" + " " + outDir;
            strcpy(mainCmd, mainTmp.c_str());
            cmdExec(mainCmd, 0);
            
            mainTmp = "rm -rf " + GSVTmp.runID + "_GSVMiningRes";
            strcpy(mainCmd, mainTmp.c_str());
            cmdExec(mainCmd, 0);
        }
        
        // pipeline ends, close log file
        
        GSVLogOut.close();
    }

    cout << "\n";
    cout << ">> " << currentDateTime().c_str() << "\tGSVMining Analysis Done!\n";

    // free memory

    delete[] mainCmd;
    
    return 0;
}

