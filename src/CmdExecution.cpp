/*
 
 CmdExecution.cpp
 
 This program is to do:
 1. run cmd
 
 Copyright (c) 2014 DeepSeq University of Nottingham. All rights reserved.
 
 Fei Sang
 DeepSeq, Queen's Medical School
 School of Life Sciences, University of Nottingham
 Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
 Email: fei.sang@nottingham.ac.uk
 
*/

#include "CurrentTime.h"
#include "CmdExecution.h"

#include <iostream>

using namespace std;

void cmdExec(char *varCmd, int varShow)
{
    if(varShow == 1)
    {
        cout << ">> " << currentDateTime().c_str() << " <<\t" << varCmd << "\n";
    }
    else
    {
        int cmdRes;
    
        cmdRes = system(varCmd);
    
        if(WEXITSTATUS(cmdRes) != 0)
        {
            cout << "ERROR!!! CMD RUN FAILED!\n";
            cout << "FAILED CMD: " << varCmd << "\n";
            exit(1);
        }
    }
}

string cmdGetOut(string varCmd, int varShow)
{
    if(varShow == 1)
    {
        cout << ">> " << currentDateTime().c_str() << " <<\t" << varCmd << "\n";
        
        string data = "test";
        
        return data;
    }
    else
    {
        string data;
        FILE * stream;
        const int max_buffer = 1024;
        char buffer[max_buffer];
    
        stream = popen(varCmd.c_str(), "r");
        if(stream)
        {
            while(!feof(stream))
            {
                if(fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
            }
            pclose(stream);
        }
    
        return data;
    }
}

