/*
 
 CurrentTime.cpp
 
 This program is to do:
 1. get current time
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
*/

#include "CurrentTime.h"

#include <iostream>

using namespace std;

string currentDateTime()
{
    time_t nowTime = time(0);
    struct tm timeStruct;
    char timeBuff[100];
    timeStruct = *localtime(&nowTime);
    strftime(timeBuff, sizeof(timeBuff), "%d-%m-%Y-%X", &timeStruct);
    
    return timeBuff;
}