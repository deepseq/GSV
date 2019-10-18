/*
 
 StringOperation.h
 
 This program is to do:
 1. Operate String
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
 including:
            split           // split string by delim(such as " .,/")
            int2str         // convert int to string
            float2str       // convert double to string
            str2int         // convert string to int
            str2float       // convert string to double
            strtoupper      // all to upper case
            strtolower      // all to lower case
 
*/

#ifndef __GSVMining__StringOperation__
#define __GSVMining__StringOperation__

#include <iostream>

#include <string>
#include <cstring>
#include <vector>
#include <sstream>
#include <algorithm>

using namespace std;

vector<string> split(string& sOStr, const char *sODel);
string int2str(int sOn);
string float2str(double sOf);
string longint2str(long int sOn);
string joinstring(vector<string> sOVec, const char *sODel);
template<class Type> string tostring(Type sOType);
int str2int(string& sOStr);
double str2float(string& sOStr);
long int str2longint(string& sOStr);
template<class ToType,class FromType> ToType strconvert(FromType sOType);
string& strtoupper(string& sOStr);
string strtoupper(string sOStr);
string& strtolower(string& sOStr);
string strtolower(string sOStr);

#endif /* defined(__GSVMining__StringOperation__) */
