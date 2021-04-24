#include <stdio.h>
#include <iostream>
#include <string.h>
#include <string>
#include <map>
#include <math.h>
#include <stdlib.h>
#include <cstdio>
#include <vector>
#include <fstream>
#include <chrono>
#include <sys/resource.h>

#include "kernel/intersection_partition.h"
#include "phasing/Phasing.h"

#define PROGRAM_BIN "main"

static const char *STRIDE_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " <command> [options]\n"  
"               phasing      test\n"
"\n";

int main(int argc, char** argv)
{
    int who = RUSAGE_SELF;
    struct rusage usage;
    int ret;
    auto begin = std::chrono::high_resolution_clock::now();

   
    if(argc <= 1)
    {
        std::cout << STRIDE_USAGE_MESSAGE;
        return 0;
    }
    
    std::string command(argv[1]);
    
    if(command=="phasing")
    {
        PhasingMain(argc - 1, argv + 1);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto welapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    cout << setiosflags(ios::left)<< setw(40) << "Time cost: " <<setprecision(2) << welapsed.count() * 1e-9 <<"s"<<endl;
    ret = getrusage(who, &usage);
    double m = (double)usage.ru_maxrss/(double)(1024*1024);
    printf("Memory usage: %lf GB\n",m);
    return 0;
}


