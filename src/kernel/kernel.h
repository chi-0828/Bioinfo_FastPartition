#ifndef KERNEL_H
#define KERNEL_H
#include <omp.h>
#include <thread>
#include <bitset>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <string>
#include <map>
#include <math.h>
#include <stdlib.h>
#include <cstdio>
#include <vector>
#include <fstream>
#include <bits/stdc++.h> 
#include <sstream>
#include <chrono>
#include "../phasing/Util.h"
#include "block.h"

using namespace std;

#define leading_bit_mask 0X7FFF
#define NEXT 1
#define PREVIOUS 0

typedef  std::map<int, vector<int>> bipartition;

struct ReadInfo{
public:
    // read_index (read ID from whatshap)
    int read_index;
    // genotype (REF or ALT)
    char genotype;
    // mapping quality
    int cost;
};

class Kernel{
public:
    Kernel(){
        positions = new vector<unsigned int>();
        snp_iterator = 0;
        inverse_state = false;
        block_num = 0;
    }
    ~Kernel(){
        delete positions;
    }
    int block_num ;
    int block_size ;
    Blockset blockset;

    // store end position of each read
    std::vector<unsigned int> read_last_pos , read_first_pos;

    // use to iterator all SNP 
    unsigned int snp_iterator; 

    //store partition result
    unordered_map <int, bipartition> current_partition,previous_partition;  

    //store reads which appear in two consecutive snp
    vector<int> forward_intersect_read,backward_intersect_read;    
    
    // store movement cost
    unordered_map <int, unsigned long> current_snp_movement,previous_snp_movement;
    // partition id -> read partition encode result
    unordered_map <int, uint64_t> encode_on_cur_partition,encode_on_pre_partition;

    // store two group of read
    unordered_map <int, uint8_t> readset; 

    // best answer for trace back
    int best_pos;

    // snp position info 
    vector<unsigned int>* positions;

    // change two tool
    //unordered_map <unsigned int , vector <pair< int, int>> >partition_move_two;

    // snp length
    unsigned int snp_size ;

    // store reads ID on each snp 
    unordered_map < unsigned int, vector<ReadInfo>> reads_on_SNP;

    //std::vector< string> all_snp;

    // be used in change_bits()
    uint64_t size_mask , partition_mask ;

    // for trace back and final partition 
    bool inverse_state ;

    // record best parition
    std::vector<unordered_map <int,int>> trace_back;
    vector<unordered_map<int,int>> trace_back_partition;
    vector<int> best_answer;
    unordered_map  <int , pair<int,int>> partition_read;

    // get read info on snp
    void get_read_info(int snp_iterator);
    
    // patition on each SNP
    void partition();

    // sorting
    void insertionSort(vector<int> &data);    
    
    // pass data to next SNP
    void pass_data();
       
    // find intersection SNP(0-1)
    void find_same_reads();
    // find intersection SNP(1-M)
    void find_same_reads_in_accumlate();
   
    // combine two SNP
    void assemble();
    void assemble_read(const int pro,int pre);

    // find the match partition at last SNP
    int find_match_encode(int currentID , uint64_t encode);
    // if no match partition , find the closest one
    void change_bits(int pos, const uint64_t &key_want_to_add);
    // get mask for XOR result
    void get_masks(int pre , int flag = 0);
    // retrun numbers of 1 in binary-n
    int count_bit(uint64_t &n);

     // make encode depend on intersection
    void build_table(int flag);
    // store the encode
    void make_encode(const int &pos,const uint64_t &enco);

    void supperread_th(int paritionID);
    //void lifehack();

    // trace back 
    // trace back -> find best path
    void trace();
    void using_best_answer_separate_two_group();
    void determine_value_depend_on_group(int &value_a , int &value_b ,int snp_id);
    void re_partition();
    bool inverse_value(int pos);
    // check if read is moved in kth parition on SNP(iterator)
    bool find_read_in_parition(int k , int iterator);

    bool __check_group(int snp_id , int read_id, int type);
};
#endif