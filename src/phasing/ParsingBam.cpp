#include "ParsingBam.h"
#include <thread>
#include <string.h>

VariantParser::VariantParser(std::string variantFile):variantFile(variantFile){
    //http://wresch.github.io/2014/11/18/process-vcf-file-with-htslib.html
    
    std::cout<< "read vcf file... \n";
    
    // counters
    int nseq = 0;
    htsFile * inf = bcf_open(variantFile.c_str(), "r");
    // read header
    bcf_hdr_t *hdr = bcf_hdr_read(inf);
    
    // report names of all the sequences in the VCF file
    const char **seqnames = NULL;
    // chromosome idx and name
    seqnames = bcf_hdr_seqnames(hdr, &nseq);
    /*
    fprintf(stderr, "Sequence names:\n");
    for (int i = 0; i < nseq; i++) {
        // bcf_hdr_id2name is another way to get the name of a sequence
        fprintf(stderr, "  [%2i] %s (bcf_hdr_id2name -> %s)\n", i, seqnames[i],bcf_hdr_id2name(hdr, i));
    }
    */
    
    std::string allSmples = "-";
    // limit the VCF data to the sample name passed in
    bcf_hdr_set_samples(hdr, allSmples.c_str(), 0);
    if (bcf_hdr_nsamples(hdr) != 1) {
        fprintf(stderr, "ERROR: please limit to a single sample\n");
    }

    // struc for storing each record
    bcf1_t *rec = bcf_init();
    
    int n    = 0;  // total number of records in file
    int nsnp = 0;  // number of SNP records in file
    int nhq  = 0;  // number of SNPs for the single sample passing filters

    // genotype data for each call
    // genotype arrays are twice as large as
    // the other arrays as there are two values for each sample
    int ngt_arr = 0;
    //int ngt     = 0;
    int *gt     = NULL;
       
    while (bcf_read(inf, hdr, rec) == 0) {
        
        n++;
        
        // this function directly looking Ref and Alt string length
        if (bcf_is_snp(rec)) {
            nsnp++;
        } else {
            continue;
        }

        //ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
        bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
        
        if ( gt[0] == 2 && gt[1] == 4) {
            // get chromosome string
            std::string chr = seqnames[rec->rid];
            // position is 0-base
            int variantPos = rec->pos;
            // get r alleles
            RefAlt tmp;
            tmp.Ref = rec->d.allele[0]; 
            tmp.Alt = rec->d.allele[1];
            // record 
            chrVariant[chr][variantPos] = tmp;

            nhq++;
            //printf("%s\t%i\t%s\t%s\n", seqnames[rec->rid], rec->pos, rec->d.allele[0], rec->d.allele[1]);
            
            //std::cout<< gt[0] << gt[1] << "\n";
            }
    }
    
    /*
    for(std::map<std::string, std::map< int, RefAlt > >::iterator chrIter = chrVariant.begin(); chrIter != chrVariant.end() ; chrIter++ ){
        //std::cout<< (*chrIter).first << "\t" << (*chrIter).second.size() << "\n";
        for(std::map< int, RefAlt >::iterator variantIter = (*chrIter).second.begin(); variantIter != (*chrIter).second.end() ; variantIter++ ){
            std::cout<< (*chrIter).first << "\t" << (*variantIter).first << "\t" << (*variantIter).second.Ref << "\t" << (*variantIter).second.Alt << "\n";
        }
    }*/
    //std::cout<< "SNP: "<< nhq << "/" << nsnp << "\n";
    
}
VariantParser::~VariantParser(){
    
}

std::map<int, RefAlt> VariantParser::getVariants(std::string chrName){
    std::map<int, RefAlt> targetVariants;
    std::map<std::string, std::map<int, RefAlt> >::iterator chrIter = chrVariant.find(chrName);
    
    if( chrIter != chrVariant.end() )
        targetVariants = (*chrIter).second;
    
    return targetVariants;
}



BamParser::BamParser(std::string inputBamFile):BamFile(inputBamFile){

}

BamParser::~BamParser(){
    
}
typedef struct samview_settings {
    void* bed;
} samview_settings_t;

void BamParser::usable_alignments(){
    
    auto region_begin = std::chrono::high_resolution_clock::now();
    // init data structure and get  core n
    htsThreadPool threadPool = {NULL, 0};
    unsigned num_cpus = std::thread::hardware_concurrency();
    
    samFile *fp_in = hts_open(BamFile.c_str(),"r"); //open bam file
	bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
	bam1_t *aln = bam_init1(); //initialize an alignment
    
    // creat thread pool
    if (!(threadPool.pool = hts_tpool_init(num_cpus))) {
        fprintf(stderr, "Error creating thread pool\n");
    }
    hts_set_opt(fp_in, HTS_OPT_THREAD_POOL, &threadPool);
    std::cout<<"Read bam file ...\n";
    ;
    char *arg[num_cpus] = {};
    long int chunk_size = 259239848/num_cpus;
    for(int i=0;i<num_cpus;i++){
        arg[i] = new char[1000];
        int start = chunk_size *i ;
        int end = chunk_size * (i+1)-1  ; // end?
        sprintf(arg[i],"%d:%d-%d",1,start ,end );
    }
    
    auto idx = sam_index_load(fp_in, (BamFile.c_str()));
    auto *iter =sam_itr_regarray(idx,bamHdr , arg , num_cpus);
    int result ;
    int count = 0;
    while ((result = sam_itr_multi_next(fp_in, iter, aln)) >= 0) {
        int flag = aln->core.flag;

        if ( aln->core.qual < 20 || //mapping quality
             (flag & 0x4)   != 0 || // read unmapped
             (flag & 0x100) != 0 || // not primary alignment
             (flag & 0x400) != 0 || // duplicate 
             (flag & 0x800) != 0){  // supplementary alignment
            continue;
        }
        
        uint8_t *q = bam_get_seq(aln); //quality string
        uint32_t len = aln->core.l_qseq; //length of the read.
        char *qseq = (char *)malloc(len);
        for(unsigned int i=0; i< len ; i++){
			qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
		}
        
        Alignment tmp;

        tmp.chr = bamHdr->target_name[aln->core.tid];
        tmp.refStart = aln->core.pos;
        tmp.qname = bam_get_qname(aln);
        tmp.qlen = aln->core.l_qseq;
        tmp.qseq = qseq;
        tmp.cigar_len = aln->core.n_cigar;
        tmp.op = (int*)malloc(tmp.cigar_len*sizeof(int));
        tmp.ol = (int*)malloc(tmp.cigar_len*sizeof(int));

        uint32_t *cigar = bam_get_cigar(aln);
        for(unsigned int k =0 ; k < aln->core.n_cigar ; k++){
            tmp.op[k] = bam_cigar_op(cigar[k]);
            tmp.ol[k] = bam_cigar_oplen(cigar[k]);
        }

        alnVec.push_back(tmp);
        count ++ ;
    }
    //std::cout<<"size"<<alnVec.size()<<std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - region_begin);
    //std::cout << "Read bam ..." <<  elapsed.count() * 1e-9 <<"s"<<std::endl;
    
   /*

	while(sam_read1(fp_in,bamHdr,aln) > 0){

        int flag = aln->core.flag;

        if ( aln->core.qual < 20 || //mapping quality
             (flag & 0x4)   != 0 || // read unmapped
             (flag & 0x100) != 0 || // not primary alignment
             (flag & 0x400) != 0 || // duplicate 
             (flag & 0x800) != 0){  // supplementary alignment
            continue;
        }
        
        uint8_t *q = bam_get_seq(aln); //quality string
        uint32_t len = aln->core.l_qseq; //length of the read.
        char *qseq = (char *)malloc(len);
        for(unsigned int i=0; i< len ; i++){
			qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
		}
        
        Alignment tmp;

        tmp.chr = bamHdr->target_name[aln->core.tid];
        tmp.refStart = aln->core.pos;
        tmp.qname = bam_get_qname(aln);
        tmp.qlen = aln->core.l_qseq;
        tmp.qseq = qseq;
        tmp.cigar_len = aln->core.n_cigar;
        tmp.op = (int*)malloc(tmp.cigar_len*sizeof(int));
        tmp.ol = (int*)malloc(tmp.cigar_len*sizeof(int));

        uint32_t *cigar = bam_get_cigar(aln);
        for(unsigned int k =0 ; k < aln->core.n_cigar ; k++){
            tmp.op[k] = bam_cigar_op(cigar[k]);
            tmp.ol[k] = bam_cigar_oplen(cigar[k]);
        }

        alnVec.push_back(tmp);

        /*
        for(std::vector<alignment>::iterator alnIter = alnVec.begin() ; alnIter != alnVec.end() ; alnIter++ ){
            std::cout<< (*alnIter).chr << "\t"
                     << (*alnIter).qname << "\t"
                     << (*alnIter).refStart << "\t"
                     << (*alnIter).cigar_len << "\t"
                     << (*alnIter).qlen << "\t";
                     //<< (*alnIter).qseq << "\t";
            for(unsigned int i = 0;i<(*alnIter).cigar_len;i++){
                std::cout<< "(" << (*alnIter).op[i] << "," << (*alnIter).ol[i] << ")";
            }
            std::cout << "\n";
        }

        std::cout << "\n";
        */        

        /*
        // debug 
        int32_t pos = aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
		char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
		//uint32_t len = aln->core.l_qseq; //length of the read.
		
		//uint8_t *q = bam_get_seq(aln); //quality string
		uint32_t q2 = aln->core.qual ; //mapping quality
		//uint32_t *cigar = bam_get_cigar(aln);
		char* qname = bam_get_qname(aln);
		//char *qseq = (char *)malloc(len);
        
        for(unsigned int i=0; i< len ; i++){
			qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
		}
		
		//printf("%s\t%d\t%s\t%d\t%d\t%d\t%s\n",chr,pos,qname,q2,flag,len,qseq);
        printf("%s\t%d\t%s\t%d\t%d\t%d\n",chr,pos,qname,q2,flag,len);
        std::cout<< "aln->core.n_cigar: " << aln->core.n_cigar << "\n";
        for(unsigned int k =0 ; k < aln->core.n_cigar ; k++){
            int op = bam_cigar_op(cigar[k]);
            int ol = bam_cigar_oplen(cigar[k]);
            
            printf("(%d,%d)",op,ol);
        }
        printf("\n");
        
	}

    std::cout<<"size "<<alnVec.size()<<std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - region_begin);
    std::cout << "Read file ori: " <<  elapsed.count() * 1e-9 <<"s"<<std::endl;**/
	bam_destroy1(aln);
	sam_close(fp_in);
}


std::vector<ReadVariant> BamParser::detect_alleles(VariantParser snpMap){
    
    std::vector<ReadVariant> reulstReadVariant;
    
    std::string currentChr;
    // use to find current chromosome
    if( alnVec.size() > 0 ){
        std::vector<Alignment>::iterator tmpIter = alnVec.begin();
        currentChr = (*tmpIter).chr;
        //std::cout<< "current chr: "<< currentChr << "\n" ;
    }
    else{
        std::cout<< "no alignment info, please check bam file.\n";
        exit(1);
    }
    
    // use chromosome to find recorded snp map
    std::map<int, RefAlt> currentVariants = snpMap.getVariants(currentChr);
    // set skip variant start iterator
    std::map<int, RefAlt>::iterator firstVariantIter = currentVariants.begin();
    if( firstVariantIter == currentVariants.end() ){
        std::cout<< "error chromosome name or empty map.\n";
        exit(1);
    }
    else{
        //std::cout<< "currentVariants size: " << currentVariants.size() << "\n";
    }    
    
    for(std::vector<Alignment>::iterator alnIter = alnVec.begin() ; alnIter != alnVec.end() ; alnIter++ ){
        ReadVariant *tmpReadResult = new ReadVariant();
        (*tmpReadResult).read_name = (*alnIter).qname;
        (*tmpReadResult).source_id = (*alnIter).chr;
        (*tmpReadResult).reference_start = (*alnIter).refStart;
        
        // Skip variants that are to the left of this read
        while( firstVariantIter != currentVariants.end() && (*firstVariantIter).first < (*alnIter).refStart )
            firstVariantIter++;
        
        //std::cout<< "first snp pos: " << (*firstVariantIter).first << "\n";
        
        // position relative to reference
        int ref_pos = (*alnIter).refStart;
        // position relative to read
        int query_pos = 0;
        // translation char* to string;
        std::string qseq = (*alnIter).qseq;
        
        // set variant start for current alignment
        std::map<int, RefAlt>::iterator currentVariantIter = firstVariantIter;
        
        // reading cigar to detect snp on this read
        for(int i = 0; i < (*alnIter).cigar_len ; i++ ){
            int cigar_op = (*alnIter).op[i];
            int length   = (*alnIter).ol[i];
            
            //std::cout<< ref_pos << "\t" << query_pos << "\t" << cigar_op << "\t" << length << "\n";
            
            // iterator next variant
            while( currentVariantIter != currentVariants.end() && (*currentVariantIter).first < ref_pos )
                currentVariantIter++;
                
            // CIGAR operators: MIDNSHP=X correspond 012345678
            // 0: alignment match (can be a sequence match or mismatch)
            // 7: sequence match
            // 8: sequence mismatch
            if( cigar_op == 0 || cigar_op == 7 || cigar_op == 8 ){
                
                while( currentVariantIter != currentVariants.end() && (*currentVariantIter).first < ref_pos + length){
                    
                    int refAlleleLen = (*currentVariantIter).second.Ref.length();
                    int altAlleleLen = (*currentVariantIter).second.Alt.length();
                    int offset = (*currentVariantIter).first - ref_pos;
                    std::string base = qseq.substr(query_pos + offset, 1);
                    int allele = -1;
                    
                    //std::cout<< "ref_pos: " << ref_pos << "\n";
                    //std::cout<< "offset: " << offset << "\n";
                    //std::cout<< "snp: " << (*currentVariantIter).second.Ref << " " << (*currentVariantIter).second.Alt << "\n";
                    //std::cout<< "base: " << base << "\n";
                    
                    if( refAlleleLen == 1 && altAlleleLen == 1){
                        if( base == (*currentVariantIter).second.Ref )
                            allele = 0;
                        else if( base == (*currentVariantIter).second.Alt )
                            allele = 1;
                    }
                    else{
                        std::cout<< "not snp: " <<(*currentVariantIter).first + 1 << "\n";
                        exit(1);
                    }
                    
                    if( allele != -1 ){
                        // record snp result
                        Variant *tmpVariant = new Variant((*currentVariantIter).first, allele, 30);
                        (*tmpReadResult).variantVec.push_back( (*tmpVariant) );
                    }
                    
                    currentVariantIter++;
                }
                query_pos += length;
                ref_pos += length;
            }
            // 1: insertion to the reference
            else if( cigar_op == 1 ){
                query_pos += length;
            }
            // 2: deletion from the reference
            else if( cigar_op == 2 ){
                ref_pos += length;
            }
            // 3: skipped region from the reference
            else if( cigar_op == 3 ){
                ref_pos += length;
            }
            // 4: soft clipping (clipped sequences present in SEQ)
            else if( cigar_op == 4 ){
                query_pos += length;
            }
            // 5: hard clipping (clipped sequences NOT present in SEQ)
            // 6: padding (silent deletion from padded reference)
            else if( cigar_op == 5 || cigar_op == 6 ){
                // do nothing
            }
            else{
                std::cout<< "alignment find unsupported CIGAR operation from read: " << (*alnIter).qname << "\n";
                exit(1);
            }
        }
        
        /*
        std::cout<< (*alnIter).chr       << "\t"
                 << (*alnIter).qname     << "\t"
                 << (*alnIter).refStart  << "\t"
                 << (*alnIter).cigar_len << "\t"
                 << (*alnIter).qlen      << "\t"
                 << (*alnIter).qseq      << "\n";
                 
        for(unsigned int i = 0;i<(*alnIter).cigar_len;i++)
            std::cout<< "(" << (*alnIter).op[i] << "," << (*alnIter).ol[i] << ")";
 
        std::cout << "\n";
        */
        
        if( (*tmpReadResult).variantVec.size() > 0 )
            reulstReadVariant.push_back((*tmpReadResult));
    }
    return reulstReadVariant;
}

