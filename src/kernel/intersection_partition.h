#ifndef INTERSECTION_PARTITION_H
#define INTERSECTION_PARTITION_H

#include <array>
#include <vector>
#include <memory>
#include <stdint.h>
#include "kernel.h"
#include "../phasing/Util.h"

class Intersection_partition{
private:
	vector<ReadVariant> * read_set;

	// following 2 data structure is order to help output result
	map<int , pair< uint8_t  , uint8_t >> result ;
	map<int , int> block_start;

    // start our phasing algorithm.
	Kernel *kernel_table;

	void compute_table();

	void writingResult(std::string vcfFile);

public:
	Intersection_partition(vector<ReadVariant> * read_set , std::vector<unsigned int> read_last_pos ,std::vector<unsigned int> read_first_pos);
 
	~Intersection_partition();

	void get_outputread(std::string vcfFile);

};

#endif
