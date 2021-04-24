#include "PhasingProcess.h"
#include "ParsingBam.h"
#include <iostream>    
#include <bits/stdc++.h>
#include <vector>
#include <map>
#include <set>
#include"../select/Readset.h"

PhasingProcess::PhasingProcess(PhasingParameters params)
{
    //std::cout << "start: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\n";
    
    VariantParser snpFile(params.snpFile);
    
    //std::cout << "parsing variant: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\n";
    //begin_time = std::clock();
    
    // TODO: each chromosome need different bam parser
    
    BamParser bamParser(params.bamFile);
    //std::cout << "create bamParser: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\n";
    //begin_time = std::clock();
    bamParser.usable_alignments();
    
    //std::cout << "usable_alignments: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\n";
    //begin_time = std::clock();
    readVariant = bamParser.detect_alleles(snpFile);
    //std::cout << "detect_alleles: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\n";

    // select 
    cout<<"select read by max-coverage 15X ... \n";
    auto begin = std::chrono::high_resolution_clock::now();
    int count = 0;
    std::vector<ReadVariant> vec_2 ;
    for(std::vector<ReadVariant>::iterator readIter = readVariant.begin() ; readIter != readVariant.end() ; ){
        if(readIter->variantVec.size() >= 2){
            count++;
            //cout<<endl;
            vec_2.push_back(*readIter);
        }
        ++readIter;
    }

    set<int>preferred_source_ids;
    Readset* readset = new Readset(vec_2);
    set<int>selected_read_set = readset->readselsect(readset,15, preferred_source_ids,true );

    
    std::vector<ReadVariant> select_reads;
    int ID = 0;
    for( std::vector<ReadVariant>::iterator readIter = readset->read.begin() ; readIter != readset->read.end() ;++readIter){

        auto target_read = selected_read_set.find((*readIter).read);

        if(target_read != selected_read_set.end()){
            read_last_pos.push_back((*readIter).variantVec.rbegin()->position );
            read_first_pos.push_back((*readIter).variantVec.begin()->position );
            //cout<<(*readIter).variantVec.end()->position<<"-"<<(*readIter).variantVec.begin()->position<<"\n";
            (*readIter).ID = ID;
            ID ++ ;
            select_reads.push_back(*readIter);
            selected_read_set.erase((*readIter).read);
        }
    }
    
    readVariant = select_reads;

    delete readset;
    readset = nullptr;
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    cout << setiosflags(ios::left)<< setw(40) << "select read : " <<setprecision(4) << elapsed.count() * 1e-9 <<"s"<<endl;
    return;
    
};

PhasingProcess::~PhasingProcess(){

};

