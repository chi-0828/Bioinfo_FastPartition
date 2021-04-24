#include <iostream>    
#include <bits/stdc++.h>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <unordered_set>
#include "Readset.h"
#include "Priorityqueue.h"
#include <algorithm>
using namespace std;
priority_type_ptr _compute_score_for_read(Readset* &readset, int index,map<int32_t,int32_t> &vcf_indices);


Priorityqueue* _construct_priorityqueue( Readset* &readset,set<int> &read_indices,map<int32_t,int32_t>  &vcf_indices){
    Priorityqueue *pqueue = new Priorityqueue();
    
    for(set<int>::iterator it =read_indices.begin();it!=read_indices.end();++it){
        //cerr<<"www\n";
        priority_type_ptr computed_score = _compute_score_for_read(readset,*it, vcf_indices);
        //cerr<<*it<<endl;
        pqueue->c_push(computed_score,*it);
        //std::sort(pqueue.begin(), pqueue.end(), mycompare);
    }
    
    return pqueue;
    
}