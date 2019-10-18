/*
 
 MappingReads.h
 
 This program is to do:
 1. map reads to reference
 2. filter mapping results
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
*/

#ifndef __GSVMining__MappingReads__
#define __GSVMining__MappingReads__

#include <iostream>

using namespace std;

void mapReads(string mRRef, string mRDir, string mRID, int mRThread);

#endif /* defined(__GSVMining__MappingReads__) */
