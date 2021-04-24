#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <vector>
#include <map>



// use for parsing
struct Variant{
    Variant(int position, int allele, int quality):position(position), allele(allele), quality(quality){};
    
    int position;
    int allele;
    int quality;
};

struct ReadVariant{
    // init function
    ReadVariant(): read_name(""), 
            mapping_quality(0), 
            source_id(""), 
            sample_id(""), 
            reference_start(0), 
            BX_tag(""){}
            
    std::string read_name;
    int mapping_quality;
    std::string source_id;
    std::string sample_id;
    int reference_start;
    std::string BX_tag;
    int read;
    int ID ;
    
    std::vector<Variant> variantVec;
};

typedef struct ReadVariant Read;

std::string getTargetString(std::string line, std::string start_sign, std::string end_sign);





#endif