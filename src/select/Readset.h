#include <iostream> 
#include <bits/stdc++.h>
#include <vector>
#include <set>
#include <map>
#include "../phasing/Util.h"
#ifndef _READSET_
#define _READSET_
using namespace std;

class Readset{
    private:
         
    public:
        Readset(vector<Read> &read){
            this->read = read;
        }
   
    map<int , int> get_all_snp(vector<Read> &read);
    void make_read_encode(vector<Read> &read,set<int> &undecided_reads);
    set<int> get_positions(map<int , int> &all_snp_records);
    
    Read get(int index,vector<Read> &read);
    set<int> readselsect(Readset* &readset, int max_cov, set<int> &preferred_source_ids, bool bridging);
    vector<Read> read;
    ~Readset(){  
        ;
    }
};
#endif     