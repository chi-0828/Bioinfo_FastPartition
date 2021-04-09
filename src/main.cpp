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

#include "../arg/argparse.hpp"
#include "intersection_partition.h"

#define PRINT(x) cout<<x<<endl

int main(int argc , char**argv){
    cxxopts::Options options("quickphasing", "Itersection pahsing");
    options.add_options()
    ("r,reference", "File name", cxxopts::value<std::string>())
    ("b,bam", "File name", cxxopts::value<std::string>())
    ("v,vcf", "File name", cxxopts::value<std::string>())
    ("0,output", "File name", cxxopts::value<std::string>())
    ("h,help", "Print usage")
    ;

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }
    string reference , bam ,vcf ,outputfile;
    if (result.count("reference"))
        reference = result["reference"].as<string>();
    else
        PRINT("No reference file");
    if (result.count("bam"))
        bam = result["bam"].as<string>();
    else{
        PRINT("Please input bam file");
        exit(1);
    }
    if (result.count("vcf"))
        vcf = result["vcf"].as<string>();
    else{
        PRINT("Please input vcf file");
        exit(1);
    }
    if (result.count("output"))
        outputfile = result["output"].as<string>();
    else{
        PRINT("Please input output file");
        exit(1);
    }

    cout<<reference<<bam<<vcf<<endl;
    //start()
}