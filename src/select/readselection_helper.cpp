#include <iostream>    
#include <bits/stdc++.h>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <unordered_set>
#include "Readset.h"
#include <algorithm>
#include "Connectcomponent.h"
#include "Priorityqueue.h"
#include <ctime> // clock 函數所需之標頭檔
using namespace std;
Priorityqueue* _construct_priorityqueue( Readset* &readset,set<int> &read_indices,map<int32_t,int32_t>  &vcf_indices);
void _slice_read_selection( Priorityqueue* &pq, map<int , int> &coverages,int  max_cov,Readset* &readset,map<int32_t,int32_t> &vcf_indices,map<int,set<int>> &variant_to_reads_map,set<int> &reads_in_slice,set<int> &reads_violating_coverage);
set<int> readselection_helper(map<int , int> &coverages, int  max_cov,Readset* &readset,map<int32_t,int32_t> &vcf_indices,map<int,set<int>> &variant_to_reads_map,set<int> &selected_reads,set<int> &undecided_reads,set<int> &positions,bool &bridging){
    Priorityqueue *pq ;
    Read read;
    //int loop = 0;
    
    while(undecided_reads.size()>0){
       
        set<int> reads_in_slice, reads_violating_coverage;
        //cerr<<"go in"<<endl;
        pq = _construct_priorityqueue(readset, undecided_reads, vcf_indices);
       // cerr<<"construct_priorityqueue done\n";
        /*for(vector<pair<vector<int>,int>>::iterator it = pq.begin();it!=pq.end();it++){
            cout<<it->second<<' '<<it->first.at(0)<<endl;
        }*/
        //cout<<endl;
        /*for(vector<pair<vector<int>,int>>::iterator it = pq.begin();it!=pq.end() ;++it){
        cout<<it->first.at(0)<<' '<<it->second<<endl;
        }*/
        clock_t start, end; // 儲存時間用的變數
        start = clock(); // 計算開始時間
        _slice_read_selection(pq, coverages, max_cov, readset, vcf_indices, variant_to_reads_map,reads_in_slice, reads_violating_coverage);
        //cerr<<"_slice_read_selection done\n";
        end = clock(); // 計算結束時間
        double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        //cout << "_slice_read_selection執行時間:" << cpu_time_used <<endl;
        
        selected_reads.insert(reads_in_slice.begin(),reads_in_slice.end());
        /*cout<<"selected_reads: "<<endl;
        
        for(set<int>::iterator it = selected_reads.begin();it!=selected_reads.end();it++){
            cout<<*it<<' ';
        }
        cout<<selected_reads.size()<<endl;
        
        cout<<endl;
        cout<<"reads_in_slice: "<<endl;
        
        for(set<int>::iterator it = reads_in_slice.begin();it!=reads_in_slice.end();it++){
            cout<<*it<<' ';
        }
        cout<<endl;
        cout<<endl;
        cout<<"reads_violating_coverage: "<<endl;
        for(set<int>::iterator it = reads_violating_coverage.begin();it!=reads_violating_coverage.end();it++){
            cout<<*it<<' ';
        }
        cout<<endl;*/
        set<int> result;
        result.clear();
        set_difference(undecided_reads.begin(), undecided_reads.end(), reads_in_slice.begin(), reads_in_slice.end(), inserter(result, result.begin()));
        undecided_reads.clear();
        undecided_reads = result;
        set<int> result2;
        set_difference(undecided_reads.begin(), undecided_reads.end(), reads_violating_coverage.begin(), reads_violating_coverage.end(), inserter(result2, result2.begin()));
        undecided_reads.clear();
        undecided_reads = result2;
        //cerr << "how much undecide read after slice"<<endl;
        //cerr << undecided_reads.size()<<endl;
        /*cout<<"Undecided reads are : ";
        for(set<int>::iterator it = undecided_reads.begin();it!=undecided_reads.end();it++)
            cout<<*it<<' ';*/

        //cout<<endl;
        //Create new component finder from reads just selected
        
        //cerr<<"minus\n";
        Connectcomponent* component_finder =new Connectcomponent ();
        component_finder->set_node(component_finder->nodes,positions);
        for(set<int>::iterator it = reads_in_slice.begin();it!=reads_in_slice.end();++it){
            Read read = readset->get(*it,readset->read);
            for(size_t i=1;i< read.variantVec.size();i++){
                component_finder->merge(component_finder->nodes,read.variantVec.at(0).position, read.variantVec.at(i).position);
            }
        }
        //cerr<<"component\n";
        set<int>bridging_reads;
        if(bridging){
            pq = _construct_priorityqueue(readset, undecided_reads, vcf_indices);
           
            pair<priority_type_ptr,item_type> entry;
            while(!pq->is_empty()){
                entry = pq->c_pop(); 

                int read_index = entry.second;
                Read read = readset->get(read_index,readset->read);
                set<int>covered_blocks;

                for(size_t i=0;i< read.variantVec.size();i++)
					covered_blocks.insert(component_finder->find(component_finder->nodes,read.variantVec.at(i).position));
                int begin = vcf_indices[read.variantVec.at(0).position];
                int end = vcf_indices[read.variantVec.at(read.variantVec.size()-1).position];
                
                int currentMax=0;
                for(int i=begin;i<=end;i++){
                    if(coverages[i]>currentMax)
                        currentMax = coverages[i];
                }
                if (currentMax>= max_cov){
		            undecided_reads.erase(read_index);
                    continue;
                }
                
                if( covered_blocks.size() < 2)
                    continue;

                bridging_reads.insert(read_index);
                selected_reads.insert(read_index);
                
                for(int i=begin;i<=end;i++){
                    coverages[i]++;
                }
                undecided_reads.erase(read_index);
                for(size_t i=1;i< read.variantVec.size();i++)
                    component_finder->merge(component_finder->nodes,read.variantVec.at(0).position, read.variantVec.at(i).position);
            
            }
        }

        //cerr << "how much undecide read after bridging"<<endl;
        //cerr << undecided_reads.size()<<endl;
        //cerr<<"bridging\n";
        //cout<<endl;
        
    }
    //cout<<selected_reads.size()<<endl;
    return selected_reads;
}