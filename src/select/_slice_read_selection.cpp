#include <iostream>    
#include <bits/stdc++.h>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <unordered_set>
#include <algorithm>
#include <iterator>  
#include"Readset.h"
#include "Priorityqueue.h"
#include <ctime> // clock 函數所需之標頭檔

using namespace std;
bool compare(pair<vector<int>,int> s1, pair<vector<int>,int> s2){
   return s1.first.at(0) < s2.first.at(0);
}
priority_type_ptr _update_score_for_reads(priority_type_ptr &former_score, Readset* &readset,int32_t index, unordered_set<int> &already_covered_variants);

void _slice_read_selection(Priorityqueue* &pq,map<int , int> &coverages,int  max_cov,Readset* &readset,map<int32_t,int32_t> &vcf_indices,map<int,set<int>> &variant_to_reads_map,set<int> &reads_in_slice,set<int> &reads_violating_coverage){


//positions of variants covered by any read selected so far
set<int> already_covered_variants;

Read extracted_read ;

pair<priority_type_ptr,item_type>  entry;
int max_item;
int pos;
unordered_set<int> variants_covered_by_this_read;

clock_t start, end1; 

while(!pq->c_is_empty()){
    variants_covered_by_this_read.clear();
    entry = pq->c_pop();
    
    //cerr<<" "<<entry.second;
    //cout<<entry.second<<' '<<entry.first.at(0)<<endl;
	max_item = entry.second;
    Read extracted_read = readset->get(max_item,readset->read);
    bool covers_new_variant = false;
    //look if positions covered by this reads are already covered or not
    for(long unsigned int i=0; i<extracted_read.variantVec.size();i++){
        pos = extracted_read.variantVec.at(i).position;
        const bool is_in =( already_covered_variants.find(pos) != already_covered_variants.end());
        if(is_in)
            continue;
        else{
            covers_new_variant = true;
            //stores the positions the read coversrue;
            variants_covered_by_this_read.insert(pos);
        }
    }
    
    //only add read if it covers at least one new variant and adding it does not violate coverage constraints
    int begin = vcf_indices[extracted_read.variantVec.at(0).position];
	int end = vcf_indices[extracted_read.variantVec.back().position];

    
    int currentMax = 0;
    for(int i=begin;i<=end;i++){
        if(coverages[i]>currentMax)
            currentMax = coverages[i];
        
    }
   
    if (currentMax>= max_cov){
		reads_violating_coverage.insert(max_item);
        continue;
    }

    else if (covers_new_variant){
        
        for(int i=begin;i<=end;i++){
            coverages[i]++;
        }
		reads_in_slice.insert(max_item);
        
        set<int> reads_whose_score_has_to_be_updated ;
        //again go over the positions in the read and add them to the already_covered_variants list

		//Only the positions in the read which cover new variants are analysed
        
        for(unordered_set<int>::iterator it=variants_covered_by_this_read.begin(); it!=variants_covered_by_this_read.end(); ++it){
                //cout<<*it<<endl;
			    already_covered_variants.insert(*it);
			    reads_whose_score_has_to_be_updated.insert(variant_to_reads_map[*it].begin(),variant_to_reads_map[*it].end());
        }
        
        set<int> selected_read_set = reads_in_slice;
        set<int> decrease_set = reads_whose_score_has_to_be_updated;
        set<int> d_set;
        
        //cout<<endl;
        //set_difference(std::begin(decrease_set), std::end(decrease_set), std::begin(selected_read_set), std::end(selected_read_set),  inserter(d_set, std::end(d_set)));  
        for(set<int ,int>::iterator it = decrease_set.begin();it!=decrease_set.end() ;++it){
            //cout<<*it<<" ";
            const bool is_in =( selected_read_set.find(*it) != selected_read_set.end());
            if(is_in){
                ;
            }
            else
                d_set.insert(*it); 
        }
        //cerr<<"d_set";
        //cerr<<' '<< d_set.size();
       //set_symmetric_difference( decrease_set.begin(), decrease_set.end(), selected_read_set.begin(), selected_read_set.end(), inserter(d_set, d_set.begin()));
        /*for(set<int ,int>::iterator it = d_set.begin();it!=d_set.end() ;++it){
        cout<<*it<<' ';
        }*/
        //cout<<endl;
        priority_type_ptr old_score;
        priority_type_ptr new_score;
        for(set<int>::iterator it=d_set.begin(); it!=d_set.end(); ++it){
            old_score= pq->c_get_score_by_item(*it);
                if (old_score != NULL){
                    new_score = _update_score_for_reads(old_score,readset,*it,variants_covered_by_this_read);
                    int in = *it;
                    pq->c_change_score(in,new_score);
                }
        }
        //cout<<endl;gggg
        
    }
        
    

}
//cerr<<endl;
}