#ifndef PHASINGPROCESS_H
#define PHASINGPROCESS_H
#include "Util.h"



struct PhasingParameters
{
    int numThreads;
    std::string snpFile;
    std::string bamFile;
    std::string phasingResult;
    std::string dotFile;
};

class PhasingProcess
{
    public:
        std::vector<unsigned int> read_last_pos , read_first_pos;
        std::vector<ReadVariant> readVariant;
        PhasingProcess(PhasingParameters params);
        ~PhasingProcess();

};


#endif