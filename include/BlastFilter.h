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

#ifndef __GSVMining__BlastFilter__
#define __GSVMining__BlastFilter__

#include <iostream>

using namespace std;

void GSVBlastFilter(string GSVBFID, int GSVBFGap, int GSVBFDist);
void GSVPrimaryFilter(string GSVPFID, int GSVPFGap, int GSVPFDist, float GSVPFRPKM, int GSVPFSpan);

#endif /* defined(__GSVMining__BlastFilter__) */
