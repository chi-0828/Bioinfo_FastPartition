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
priority_type_ptr _update_score_for_reads(priority_type_ptr &former_score, Readset* &readset,int32_t index, unordered_set<int> &already_covered_variants){
    //updatest the score of the read, depending on how many reads are already covered
    int first_score = former_score->at(0);
    int second_score = former_score->at(1);
	int quality = former_score->at(2);
	Read read = readset->get(index,readset->read);
    for (int i=0;i< read.variantVec.size();i++){
        unordered_set<int>::iterator it = already_covered_variants.find(read.variantVec.at(i).position);
        if (it == already_covered_variants.end()){
            first_score--;
        }
            
    }
    priority_type_ptr result = new priority_type();
	result->push_back(first_score);
	result->push_back(second_score);
	result->push_back(quality);
	return result;
    
}