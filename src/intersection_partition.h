#ifndef INTERSECTION_PARTITION_H
#define INTERSECTION_PARTITION_H

#include <array>
#include <vector>
#include <memory>
#include "columniterator.h"
#include "entry.h"
#include "read.h"
#include "readset.h"
#include "kernel.h"

class Intersection_partition{
private:
	ReadSet* read_set;
	// stores the sample index for each read
	std::vector<unsigned int> read_sources;
	
	ColumnIterator input_column_iterator;

    // start our phasing algorithm. by GARY 2021-1-13  
    // recommend rename snp to all_snp_table or full_snp_table
	Kernel kernel_table;

	void compute_table();

	std::unique_ptr<std::vector<unsigned int> > extract_read_ids(const std::vector<const Entry *>& entries);

public:
	Intersection_partition(ReadSet* read_set, const vector<unsigned int>* positions);
 
	~Intersection_partition();

	std::pair<Read*,Read*> get_outputread();

};

#endif
