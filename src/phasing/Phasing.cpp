#include "Phasing.h"
#include "PhasingProcess.h"
#include "Util.h"
#include "../kernel/intersection_partition.h"
#include <getopt.h>


#define SUBPROGRAM "Phasing"

static const char *CORRECT_USAGE_MESSAGE =
"Usage: "  " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"      --help                     display this help and exit.\n"
"      -s  --snp-file=NAME        input SNP file.\n"
"      -b  --bam-file=NAME        input bam file.\n"
"      -d, --dot-prefix=NAME      generate dot file. \n"
"      -o, --out-prefix=NAME      phasing result.\n"
"      -t, -threads               thread. \n"
"\n";


static const char* shortopts = "s:b:d:o:t:";

enum { OPT_HELP = 1};//, OPT_ALL, OPT_IGNORE};

static const struct option longopts[] = { 
    { "help",              no_argument,       NULL, OPT_HELP },
    { "snp-file",          required_argument, NULL, 's' },
    { "bam-file",          required_argument, NULL, 'b' },
    { "dot-prefix",        no_argument, NULL, 'd' },
    { "out-prefix",        no_argument, NULL, 'o' },
    { "threads",           no_argument, NULL, 't' },
//    { "ignore-noise",      no_argument,       NULL, OPT_IGNORE},
//    { "all",               no_argument,       NULL, OPT_ALL },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static int numThreads = 1;
    static std::string snpFile="";
    static std::string bamFile="";
    static std::string phasingResult="default";
    static std::string dotFile="";
}

void PhasingOptions(int argc, char** argv)
{
    optind=1;    //reset getopt
    
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
        case 's': arg >> opt::snpFile; break;
        case 'b': arg >> opt::bamFile; break;
        case 'o': arg >> opt::phasingResult; break;
        case 'd': arg >> opt::dotFile; break;
        case OPT_HELP:
            std::cout << CORRECT_USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 0 )
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    }

    if( opt::snpFile != "")
    {
        std::ifstream openFile( opt::snpFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cout<< "File " << opt::snpFile << " not exist.\n\n";
            die = true;
        }
    }
    
    if( opt::bamFile != "")
    {
        std::ifstream openFile( opt::bamFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cout<< "File " << opt::bamFile << " not exist.\n\n";
            die = true;
        }
    }

    if (die)
    {
        std::cout << "\n" << CORRECT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }


    
}

int PhasingMain(int argc, char** argv)
{

    PhasingParameters ecParams;
    // set parameters
    PhasingOptions(argc, argv);
    // no file in command line

    ecParams.numThreads=opt::numThreads;
    ecParams.snpFile=opt::snpFile;
    ecParams.bamFile=opt::bamFile;
    ecParams.phasingResult=opt::phasingResult;
    ecParams.dotFile=opt::dotFile;
    
    PhasingProcess processor(ecParams);
    vector<ReadVariant>* readVariant = new vector<ReadVariant>(processor.readVariant.begin() ,processor.readVariant.end()) ;
    
    // start phasing algorithim
    Intersection_partition *intersection_partition = new Intersection_partition(readVariant ,processor.read_last_pos , processor.read_first_pos );
    // get final output
    intersection_partition->get_outputread(ecParams.snpFile);

    return 0;
}