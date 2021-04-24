#ifndef BLOCK_H
#define BLOCK_H
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <assert.h>
using namespace std;

class Block{
private:
    int start_pos ;
    int end_pos;
    int ID;
    int len;
public:
    Block(int start_pos , int ID ):
        start_pos(start_pos),
        ID(ID){};
    
    int getID();

    int get_start();

    int get_end();

    int get_length();

    void set_end(int end_pos);

    void set_len(int len);
    string toString();

    friend std::ostream& operator<<(std::ostream& out, const Block& e);
};
class Blockset{
private:
    vector <Block*> overall_block;
public:
    ~Blockset(){
        for(auto b : overall_block){
            delete[] b;
        }
    }
    string find_largest_block();

    int get_size();

    void add_block(Block *block , int pre_end = 0 , int len = 0);

    Block* get(int pos);

    friend std::ostream& operator<<(std::ostream& out, const Blockset& e);
};
#endif