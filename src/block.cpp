#include "block.h"
//Blockset
string Blockset::find_largest_block(){
    int ID = -1;
    int max_len =0;
    for(auto &block : overall_block){
        int len = block->get_length();
        if( len > max_len ){
            ID = block->getID();
            max_len = len;
        }
    }
    assert(ID >= 0);
    return "Largest " + overall_block[ID]->toString();
}
int Blockset::get_size(){
    return overall_block.size();
}
void Blockset::add_block(Block *block , int pre_end, int len ){
    int ID = block->getID();
    
    if(ID != -1){
        //cout<<"ID: "<<ID<<endl;
        //cout<<"get_start: "<<block->get_start()<<endl;
        overall_block.push_back(block);
        ID = this->get_size()-1;
    }else{
        ID = this->get_size();
    }

    
    if(ID != 0){
        // setting previous block end position
        overall_block[ID-1]->set_end(pre_end);
        overall_block[ID-1]->set_len(len);
    }
}

//Block::
int Block::getID(){
    return this->ID;
}

int Block::get_start(){
    return this->start_pos;
}

int Block::get_end(){
    return this->end_pos;
}

int Block::get_length(){
    return this->len;
}
    
void Block::set_end(int end_pos){
    this->end_pos = end_pos;
}
void Block::set_len(int len){
    this->len = len;
}
std::ostream& operator<<(std::ostream& out, const Blockset& b) {
    for(auto block : b.overall_block)
	    out << block->toString() <<endl ;
    
	return out;
}
string Block::toString() {
    ostringstream out ;
    out << "block(ID: "<<this->ID<<") contains "<<this->len<<" variants between position "<<this->start_pos<<" and "<<this->end_pos;
	return out.str();
}