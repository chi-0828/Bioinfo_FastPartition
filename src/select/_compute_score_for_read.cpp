#include <iostream>    
#include <bits/stdc++.h>
#include <vector>
#include <map>
#include <algorithm>
#include "Readset.h"
using namespace std;
#include "Priorityqueue.h"
priority_type_ptr _compute_score_for_read(Readset* &readset, int index,map<int32_t,int32_t> &vcf_indices){
    /*Compute the score for one read, assuming no other reads have been selected so far
	(after selecting reads, scores can be updated using _update_score_for_reads).
	We use the following scoring scheme: (new - gaps, total - bad, min(qualities)),
	where "new" is the number of variants covered by this read and no other selected read (so far),
	"gaps" is the number of variants overlapped by (physical) fragment, but not covered by the sequenced part of the read,
	"total" is the total number of variants covered by the read, and "min(qualities)" is the minimum
	over all base qualities at variant positions.
	*/
    Read read = readset->get(index,readset->read);
    int32_t min_quality = -1;
	int32_t good_score = 0;
	int32_t bad_score = 0;
	int32_t quality = -1;
	int32_t pos = -1;
    vector<int> covered_variants;
    for(size_t i=0; i<read.variantVec.size();i++){
        int quality = read.variantVec.at(i).quality;
        pos = read.variantVec.at(i).position;
        if (!i)
			min_quality = quality;
		else
			min_quality = min(min_quality, quality);
        
        int32_t variant_covered = vcf_indices[pos];
        if(variant_covered){
            covered_variants.push_back(variant_covered);
            good_score++;
        }
    }
    if(covered_variants.size()!=(covered_variants[covered_variants.size() - 1] - covered_variants.front() + 1))
        bad_score = (covered_variants.at(covered_variants.size() - 1) - covered_variants.front() + 1)-covered_variants.size();
        
    priority_type_ptr result = new priority_type();
    result->push_back(good_score - bad_score );
	result->push_back(good_score - bad_score );
	result->push_back(min_quality);
	return result;
    
    }

