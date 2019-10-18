/*
 
 CoverageFilter.h
 
 This program is to do:
 1. fitler false positive results
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
*/

#ifndef __GSVMining__CoverageFilter__
#define __GSVMining__CoverageFilter__

#include <iostream>

using namespace std;

void GSVCovFilter(string GSVCFID, int GSVCFDist, int GSVCFLen, long int GSVCFTotalChr, long int GSVCFTotalContig);

#endif /* defined(__GSVMining__CoverageFilter__) */
