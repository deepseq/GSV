/*
 
 CmdExecution.h
 
 This program is to do:
 1. run cmd
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
 */

#ifndef __GSVMining__CmdExecution__
#define __GSVMining__CmdExecution__

#include <iostream>

using namespace std;

void cmdExec(char *varCmd, int varShow);
string cmdGetOut(string varCmd, int varShow);

#endif /* defined(__GSVMining__CmdExecution__) */

