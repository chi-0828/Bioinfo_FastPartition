#include <iostream>    
#include <bits/stdc++.h>
#include <vector>
#include <set>
#include <map>
#include "Priorityqueue.h"
#include "Readset.h"

using namespace std;

void Readset::make_read_encode(vector<Read> &read,set<int> &undecided_reads ){
    int i=0;
    for(std::vector<Read>::iterator readIter = read.begin() ; readIter !=read.end() ; readIter++ ){
        undecided_reads.insert(i);
        readIter->read = i;
        i++;
    }
}

map<int , int> Readset::get_all_snp(vector<Read> &read){
    map<int , int> all_snp_records;
    for(std::vector<ReadVariant>::iterator readIter = read.begin() ; readIter != read.end() ; readIter++ )
         for(std::vector<Variant>::iterator variantIter = (*readIter).variantVec.begin(); variantIter != (*readIter).variantVec.end() ; variantIter++ )
             all_snp_records.insert(std::make_pair(variantIter->position , 0));
    
    return all_snp_records;

}

Read Readset::get(int index,vector<Read> &read){
    Read extract_read = read.at(index);
    //Read get 
    
    /*for(std::vector<Read>::iterator readIter = read.begin() ; readIter !=read.end() ; readIter++ ){
        if(readIter->read == index){
            extract_read = *readIter;
            break;
        }
    }*/
    return extract_read;
}
set<int> Readset::get_positions(map<int , int> &all_snp_records){
    
    set<int> positions;
    for(auto &snp: all_snp_records){
        positions.insert(snp.first);
    }
    return positions;

}
void  _construct_indexes(map<int , int> all_snp_records,Readset *readset,set<int> preferred_source_ids,set<int> &positions,map<int32_t,int32_t> &vcf_indices,map<int,set<int>> &variant_to_reads_map,set<int>  &preferred_reads);
set<int> readselection_helper(map<int , int> &coverages, int  max_cov,Readset* &readset,map<int32_t,int32_t> &vcf_indices,map<int,set<int>> &variant_to_reads_map,set<int> &selected_reads,set<int> &undecided_reads,set<int> &positions,bool &bridging);

set<int> Readset::readselsect(Readset* &readset, int max_cov, set<int> &preferred_source_ids, bool bridging){
clock_t start, end;
start = clock(); // 計算開始時間
set<int> positions;
map<int32_t,int32_t> vcf_indices;
map<int,set<int>> variant_to_reads_map;
set<int> preferred_reads;
set<int> undecided_reads ;
set<int> selected_reads;
map<int , int>all_snp_records = readset->get_all_snp(readset->read);
map<int , int>coverages;
readset->make_read_encode(readset->read,undecided_reads);
_construct_indexes(all_snp_records,readset,preferred_source_ids,positions,vcf_indices,variant_to_reads_map, preferred_reads);

for(long unsigned int i=0;i<undecided_reads.size();i++){
    coverages.insert(make_pair(i,0));
}
//cout<<coverages.size();

//map<int , int>coverages = all_snp_records;




for(vector<ReadVariant>::iterator it = readset->read.begin();it!=readset->read.end();++it){
    if(it->variantVec.size()<2){
        cout<<"Dectect length fail!"<<endl;
        exit(1);
    }
}
/*

if (preferred_reads.size() > 0){
    set<int> selected_preferred_reads ;
	selected_preferred_reads = readselection_helper(coverages, max_cov, readset, vcf_indices, variant_to_reads_map, selected_reads, preferred_reads, positions, bridging);
	selected_reads.insert(selected_preferred_reads);
    set<int> result;
    result.clear();
    set_symmetric_difference(undecided_reads.begin(), undecided_reads.end(), preferred_reads.begin(), preferred_reads.end(), inserter(result, result.begin()));
    undecided_reads.clear();
    undecided_reads = result;
}*/
selected_reads = readselection_helper(coverages, max_cov, readset, vcf_indices, variant_to_reads_map, selected_reads, undecided_reads, positions, bridging);
//cout<<selected_reads.size()<<endl;
/*for(std::set<int>::iterator readIter = selected_reads.begin() ; readIter !=selected_reads.end() ; readIter++ ){
        cout<<*readIter<<' '<<endl;

    }
    cerr<<endl;*/
end = clock(); // 計算結束時間
double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
//cout << "_slice_read_selection執行時間:" << cpu_time_used <<endl;

return selected_reads;


}