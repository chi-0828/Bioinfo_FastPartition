#ifndef PARSINGBAM_H
#define PARSINGBAM_H

#include "Util.h"
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/kbitset.h>
#include <htslib/vcf.h>
#include <htslib/thread_pool.h>
#include <htslib/vcfutils.h>

struct RefAlt{
    std::string Ref;
    std::string Alt;
};

class VariantParser{
    
    private:
        std::string variantFile;
        // std::map< chromosome , std::map< position , RefAlt >  >
        // position is 0-base
        std::map<std::string, std::map<int, RefAlt> > chrVariant;
    
    public:
        VariantParser(std::string inputBamFile);
        ~VariantParser();
        
        std::map<int, RefAlt> getVariants(std::string chrName);

};

struct Alignment{
    std::string chr;
    std::string qname;
    int refStart;
    int qlen;
    char *qseq;
    int cigar_len;
    int *op;
    int *ol;
};

class BamParser{
    
    private:
        std::string BamFile;
        std::vector<Alignment> alnVec;
    
    public:
        BamParser(std::string inputBamFile);
        ~BamParser();
        
        void usable_alignments();
        std::vector<ReadVariant> detect_alleles(VariantParser snpMap);
};






#endif