/*
 
 ReadsAssembly.cpp
 
 This program is to do:
 1. assemble unmapped reads
 2. summary assembly
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
*/

#ifndef __GSVMining__ReadsAssembly__
#define __GSVMining__ReadsAssembly__

#include <iostream>

using namespace std;

long int readsAss(string rALen, string rAIns, string rAsdIns, string rARef, string rAID, long int rANum, int rAGap, int rADist, int rAThread);

#endif /* defined(__GSVMining__ReadsAssembly__) */
