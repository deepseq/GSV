/*
 
 StringOperation.cpp
 
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

#include "StringOperation.h"

#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <sstream>
#include <algorithm>

using namespace std;

// sOStr string to be splited
// sODel delimiter, const char*, just like " .,/", white space, dot, comma, splash
// return a string vector saved all the splited words

vector<string> split(string& sOStr, const char *sODel)
{
    char *strOperStr, *strOperP;
    vector<string> strOperRes;
    strOperStr = new char[sOStr.size()+1];
    strcpy(strOperStr,sOStr.c_str());
    strOperP = strtok(strOperStr,sODel);
    while(strOperP!=NULL)
    {
        strOperRes.push_back(strOperP);
        strOperP = strtok(NULL,sODel);
    }
    return strOperRes;
}

// convert a integer into string through stringstream
// return the string form of n

string int2str(int sOn)
{
    stringstream strOperSS;
    string strOperStr;
    strOperSS << sOn;
    strOperSS >> strOperStr;
    return strOperStr;
}

string float2str(double sOf)
{
    stringstream strOperSS;
    string strOperStr;
    strOperSS << sOf;
    strOperSS >> strOperStr;
    return strOperStr;
}

string longint2str(long int sOn)
{
    stringstream strOperSS;
    string strOperStr;
    strOperSS << sOn;
    strOperSS >> strOperStr;
    return strOperStr;
}

// join a vector
// return the string

string joinstring(vector<string> sOVec, const char *sODel)
{
    stringstream strOperStr;
    int sOi;
    
    for(sOi = 0; sOi < int(sOVec.size()); sOi++)
    {
        if(sOi != 0){strOperStr << sODel;}
        strOperStr << sOVec[sOi];
    }
    
    return strOperStr.str();
}

// convert something to string form through stringstream
// Type Type can be int,float,double
// return the string form of param a

template<class Type> string tostring(Type sOType)
{
    stringstream strOperSS;
    string strOperStr;
    strOperSS << sOType;
    strOperSS >> strOperStr;
    return strOperStr;
}

// @brief convert string to int by atoi
// @return the integer result

int str2int(string& sOStr)
{
    return atoi(sOStr.c_str());
}

double str2float(string& sOStr)
{
    return atof(sOStr.c_str());
}

long int str2longint(string& sOStr)
{
    return atol(sOStr.c_str());
}

// do string convert through stringstream from FromType to ToType
// ToType target type
// FromType source type
// t to be converted param
// return the target form of param t

template<class ToType,class FromType> ToType strconvert(FromType sOType)
{
    stringstream strOperSS;
    ToType strOperType;
    strOperSS << sOType;
    strOperSS >> strOperType;
    return strOperType;
}

// convert string to upper case throught transform method, also can use transform method directly
// return the upper case result saved still in s

string& strtoupper(string& sOStr)
{
    transform(sOStr.begin(),sOStr.end(),sOStr.begin(),::toupper);
    return sOStr;
}

// convert string to upper case through toupper, which transform a char into upper case
// return the upper case result string

string strtoupper(string sOStr)
{
    string strOperT = sOStr;
    int sOi = -1;
    while(strOperT[sOi++])
    {
        strOperT[sOi] = toupper(strOperT[sOi]);
    }
    return strOperT;
}

// convert string to lower case throught transform method, also can use transform method directly
// return the lower case result saved still in s

string& strtolower(string& sOStr)
{
    transform(sOStr.begin(),sOStr.end(),sOStr.begin(),::tolower);
    return sOStr;
}

// convert string to lower case through tolower, which transform a char into lower case
// return the lower case result string

string strtolower(string sOStr)
{
    string strOperT = sOStr;
    int sOi = -1;
    while(strOperT[sOi++])
    {
        strOperT[sOi] = tolower(strOperT[sOi]);
    }
    return strOperT;
}
