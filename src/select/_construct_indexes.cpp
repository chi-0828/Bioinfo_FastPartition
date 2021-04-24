#include <iostream>    
#include <bits/stdc++.h>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <unordered_set>
#include "Readset.h"
#include <algorithm>

using namespace std;
void  _construct_indexes(map<int , int> all_snp_records,Readset *readset,set<int> preferred_source_ids,set<int> &positions,map<int32_t,int32_t> &vcf_indices,map<int,set<int>> &variant_to_reads_map,set<int>  &preferred_reads){
//The parameter readset: is the given ReadSet and returns, all possible variant positions, the vcf_index_ mapping and the variant_to_reads_map

positions = readset->get_positions(all_snp_records);

set<int>::iterator it;
int32_t encode = 0;

for(it=positions.begin();it!= positions.end();++it){
    set<int> l;
    variant_to_reads_map.insert( make_pair(*it,l));
    vcf_indices.insert( make_pair(*it, encode));
    encode++;
}
    for(std::vector<ReadVariant>::iterator readIter = readset->read.begin() ; readIter !=readset->read.end() ; readIter++ ){
        if (!preferred_source_ids.empty()){
            const bool is_in = preferred_source_ids.find(readIter->read) != preferred_source_ids.end();
            if(is_in){
                preferred_reads.insert(readIter->read);
            }
        }
        for(std::vector<Variant>::iterator variantIter = (*readIter).variantVec.begin(); variantIter != (*readIter).variantVec.end() ; variantIter++ ){
            variant_to_reads_map[variantIter->position].insert(readIter->read);
        }
    }
}