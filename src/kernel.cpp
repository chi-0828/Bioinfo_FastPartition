#include "kernel.h"

void Kernel::partition(){    

    //#pragma omp declare reduction (merge : std::map<int,unsigned int> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

    //#pragma omp parallel for 
    for(size_t partition_iterator=0;partition_iterator < reads_on_SNP.at(positions->at(snp_iterator)).size()+1;++partition_iterator){   
        auto snp_pos = positions->at(snp_iterator);
        // This table give the space of storing each partition for genotype 0 and 1
        bipartition table;
       
        // This vector give the space for 'table' to store each partition
        vector<int>tmp;
        table.insert(pair<int, vector<int>>(0,tmp));
        table.insert(pair<int, vector<int>>(1,tmp));

        // insert table
        current_partition.insert(pair<int,bipartition>(partition_iterator,table));
        size_t The_read_on_SNP=0;

        // to record the movement cost
        if(partition_iterator==0){
            current_snp_movement.insert(pair<int, int>(partition_iterator,0));
        }
        else{
            current_snp_movement.insert(pair<int, int>(partition_iterator,1));
        }

        // Put the read one by one to the bipartition
        while(The_read_on_SNP < reads_on_SNP[snp_pos].size()){
            // Since we want to store all sort of bipartition ,
            // so if The_read_on_SNP == partition_iterator-1 ,
            // it means we will change this read to another bipartition
            // else it will remain same in the bipartition
            // For instance , The read (4 5 6) (7 8).If The_read_on_SNP == 1
            // it will let 4(first one) change down  and the reads 5 6 7 8 remain same in bipartition
            // thus, the final bipartition will become (5 6) (4 7 8) 
            //current_partition[partition_iterator][0].reserve(reads_on_SNP[snp_pos].size());
            //current_partition[partition_iterator][1].reserve(reads_on_SNP[snp_pos].size());
            //The_read_on_SNP != partition_iterator-1 ,so remain same ,genotype 0 goes partition number 0, genotype 1 goes partition number 1
            if(The_read_on_SNP != partition_iterator-1){
                if(reads_on_SNP[snp_pos][The_read_on_SNP].genotype=='0')
                    current_partition[partition_iterator][0].push_back(reads_on_SNP[snp_pos][The_read_on_SNP].read_index);
                else
                    current_partition[partition_iterator][1].push_back(reads_on_SNP[snp_pos][The_read_on_SNP].read_index); 
            }

            //The_read_on_SNP == partition_iterator-1 ,so change ,genotype 0 goes partition number 1, genotype 1 goes partition number 0
            else{
                if(reads_on_SNP[snp_pos][The_read_on_SNP].genotype=='1')
                    current_partition[partition_iterator][0].push_back(reads_on_SNP[snp_pos][The_read_on_SNP].read_index);
                else
                    current_partition[partition_iterator][1].push_back(reads_on_SNP[snp_pos][The_read_on_SNP].read_index); 
            }
            ++The_read_on_SNP;
        }
        
    } 

}



void Kernel::assemble_read(const int pro,int pre){
    //put reads together
    
    // start assemble
    for(size_t i=0;i<previous_partition[pre][0].size();++i){
        // ignore read that has ended
        if(read_last_pos.at(previous_partition[pre][0][i]) < positions->at(snp_iterator)){
            //cout<<"ignore R:"<<previous_partition[pre][0][i]<<endl;
            continue;
        }
        if(std::find(backward_intersect_read.begin(),backward_intersect_read.end(),previous_partition[pre][0][i]) == backward_intersect_read.end()){
            if(trace_back_partition.at(trace_back_partition.size()-1).at(pro))
                current_partition[pro][1].push_back(previous_partition[pre][0][i]);

            else
                current_partition[pro][0].push_back(previous_partition[pre][0][i]);
        }
    }
    for(size_t i=0;i<previous_partition.at(pre).at(1).size();++i){
        // ignore read that has ended
        if(read_last_pos.at(previous_partition[pre][1][i]) < positions->at(snp_iterator)){
            //cout<<"ignore R:"<<previous_partition[pre][1][i]<<endl;
            continue;
        }
        if(std::find(backward_intersect_read.begin(),backward_intersect_read.end(),previous_partition[pre][1][i]) == backward_intersect_read.end()){
            if(trace_back_partition.at(trace_back_partition.size()-1).at(pro))
                current_partition[pro][0].push_back(previous_partition[pre][1][i]);

            else
                current_partition[pro][1].push_back(previous_partition[pre][1][i]);
        }
    }
    //final -> sort
    insertionSort(current_partition[pro][1]);
    insertionSort(current_partition[pro][0]);
}
int Kernel::find_match_encode(int currentID , uint64_t encode){
    int match_id = -1; 
    auto min_movement = previous_snp_movement[0];  //to find min cost
    int inverse = 0;
    for(auto const& imap: encode_on_pre_partition){
        get_masks(imap.first , 1);
        if(imap.second == (encode&size_mask) || imap.second == ((~encode)&size_mask)){
            if(previous_snp_movement[imap.first] <= min_movement){
                if(previous_snp_movement[imap.first] == min_movement){
                    // mutiple choices  
                    //mutiple_choose.at(snp_iterator).at(currentID).push_back(match_id);
                    //////////////////////////////////////////////////////////////////
                }
                match_id = imap.first;
                min_movement = previous_snp_movement[imap.first];
                if(imap.second == (encode&size_mask))
                    inverse = 0;
                else
                    inverse = 1;
            }
        }
    }
    if(match_id != -1){
        int size = trace_back.size()-1;
        trace_back[size].insert(pair<int,int>(currentID,match_id));
        if(inverse == 1){
            trace_back_partition[size].insert(pair<int,int>(currentID,1));
        }
        else{
            trace_back_partition[size].insert(pair<int,int>(currentID,0));
        }
        //cout<<"inverse : "<<inverse<<endl;
    }
    return match_id;
}
void Kernel::assemble()
{
    unordered_map<int ,int> record;
    trace_back.push_back(record);
    unordered_map<int ,int> record_partition;
    trace_back_partition.push_back(record_partition);

    for(auto const& imap: encode_on_cur_partition){
        //cout<<"cur pos"<<imap.first<<endl;
        // store info for debug (mutiple path)
        //vector <int> temp;
        //mutiple_choose.at(snp_iterator).insert(pair<int,vector<int>>(imap.first,temp));
        ///////////////////////////////////////////////////////////////////////////////
        int pre = find_match_encode(imap.first,imap.second);
        // if match successfully
        // put reads together
        if(pre != -1){
            //cout<<pre <<" match "<< imap.first<<endl;
            assemble_read(imap.first,pre);
            current_snp_movement[imap.first] += previous_snp_movement[pre];
            //cout<<std::bitset<16>(encode_on_pre_partition.at(pre)) << "\n" <<std::bitset<16>(encode_on_cur_partition.at(imap.first))<<endl;
            //cout<< previous_snp_movement.at(pre) << " : " <<current_snp_movement.at(imap.first)<<endl;
            //cout<<"pre: "<<pre<<endl;;
            //if(snp_iterator>20)
            //getchar();
        }
        // match failed , using change bits to find another encode
        else{
            //cout<<" change "<< imap.first<<endl;
            change_bits(imap.first,imap.second);
        }
    }
    bool new_block = true;
    for(auto &read : current_partition[0][0]){
        if(read_first_pos[read] != positions->at(snp_iterator))
            new_block = false;
    }
    for(auto &read : current_partition[0][1]){
        if(read_first_pos[read] != positions->at(snp_iterator))
            new_block = false;
    }
    if(new_block){
        blockset.add_block(new Block(positions->at(snp_iterator),block_num) , positions->at(snp_iterator-1) , block_size+1);
        block_num++;
        block_size = 1;
    }else{
        block_size ++ ;
    }
}
inline uint64_t POW(const unsigned int &e){
    uint64_t pow = 1 << e ; 
    return  pow;
}
void Kernel::build_table(int flag){
    encode_on_cur_partition.clear();
    vector<int> same;
    if(flag == 1){
        same = forward_intersect_read;
    }
    else{
        same = backward_intersect_read;
    }
    size_t i =0;
    //cout<<flag<<" intersection size: "<<same.size()<<endl;
    while (i < current_partition.size()){
        uint64_t encode = 0;

        //build encode
        for(size_t k=0;k<same.size();k++){
            //cout<<k<<" : "<<same.at(k)<<endl;
            if(std::find(current_partition[i][0].begin(), current_partition[i][0].end(), same.at(k)) != current_partition[i][0].end()){
                ;
            }
            else{
                encode += POW(same.size()-k-1);
            }
            //cout<<"build :encode for "<<i <<" parition :\n"<<bitset<16>(encode)<<endl;
        }
        
        
        //put in record map
        encode_on_cur_partition.insert(pair<int,uint64_t>(i,encode));
        ++i;
    }
    //getchar();
}

inline int Kernel::count_bit(uint64_t &n){
    int count = 0 ; 
    while (n) { 
        count += (n & 1); 
        n >>= 1; 
        // Force the first bit from the left to 0
        // 32767 is a number for 16 bits operation , it can be adjust for select for reads
        n = n & leading_bit_mask;
    } 
    return count;
}
inline void Kernel::get_masks(int pre , int flag){
    // a mask to avoid out of range
    int size = backward_intersect_read.size();

    size_mask = POW(size) - 1;

    // for better performance , find match encode don't need following action
    if(flag)
        return ; 
    
    if(pre == 0)
        partition_mask = 0XFFFF;
    else{
        for(int i = 0; i< size; ++i ) {
            if(backward_intersect_read[i] == reads_on_SNP[positions->at(snp_iterator-1)][pre-1].read_index){
                partition_mask = POW(size-i-1);
                partition_mask = ~partition_mask;
                break;
            }
        }
    }
    
}

void Kernel::change_bits(int pos, const uint64_t &key_want_to_add){

    auto min = 100000000;
    auto cost = encode_on_pre_partition[0];
    int pre;
    int change_value = 0;
    for(auto const& imap: encode_on_pre_partition){
        get_masks(imap.first);
        int state = 0;

        // to XOR encdoe
        uint64_t n = imap.second ^ key_want_to_add;
        
        // to make sure size range
        n = n & size_mask;
        
        // avoid moving a read more than once
        n = n & partition_mask;
        
        int diff = count_bit(n);


        // to XOR ~(encode)
        uint64_t n2 = imap.second ^ ~(key_want_to_add);
        
        // to make sure size range
        n2 = n2 & size_mask;
        ;
        // avoid moving a read more than once
        n2 = n2 & partition_mask;
        
        int diff2 = count_bit(n2);  

        // comapre 2 XOR result and choose the better one 
        if (diff2 < diff){
            diff = diff2; 
            state = 1;
        }
        if(diff <= min){
            if(diff == min){
                if(cost>=previous_snp_movement[imap.first]){
                    min = std::move(diff);
                    cost = previous_snp_movement[imap.first];
                    pre = std::move(imap.first);
                    change_value = state;
                }
            }
            else {
                min = std::move(diff);
                cost = previous_snp_movement[imap.first];
                pre = std::move(imap.first);
                change_value = state;
            }
        } 

    }
    int size = trace_back.size()-1;
    trace_back[size].insert(pair<int,int>(pos,pre));
    if(change_value)
        trace_back_partition[size].insert(pair<int,int>(pos,1));
    else
        trace_back_partition[size].insert(pair<int,int>(pos,0));
    assemble_read(pos,pre);

    //get_masks(pre);
    //cout<<"mask : \n"<<std::bitset<16>(size_mask) <<"\n"<<std::bitset<16>(partition_mask)<<endl;;
    current_snp_movement[pos] += previous_snp_movement[pre];
    current_snp_movement[pos] += min;
    /*cout<< "pre : " <<pre << " , diff : "<<min<<endl; 
    if(trace_back_partition.at(size).at(pos))
        cout<<std::bitset<16>(encode_on_pre_partition.at(pre)) << "\n" <<std::bitset<16>(~key_want_to_add)<< "\n" ;
    else
        cout<<std::bitset<16>(encode_on_pre_partition.at(pre)) << "\n" <<std::bitset<16>(key_want_to_add)<< "\n" ;
    cout<< previous_snp_movement.at(pre) << " : " <<current_snp_movement.at(pos)<<endl; 
    cout<< " inverse : " <<trace_back_partition.at(size).at(pos)<<endl; 
    //if(snp_iterator>20)
    getchar();*/
}

    
void Kernel::find_same_reads(){
    forward_intersect_read.clear();
    //find same
    unsigned int now_snp_name = positions->at(snp_iterator);
    unsigned int next_snp_name = positions->at(snp_iterator+1);
    
    size_t i =0,j=0;
    while(i <reads_on_SNP.at(now_snp_name).size() && j <reads_on_SNP.at(next_snp_name).size()) {
        if(reads_on_SNP.at(now_snp_name).at(i).read_index == reads_on_SNP.at(next_snp_name).at(j).read_index){
            forward_intersect_read.push_back(move(reads_on_SNP.at(next_snp_name).at(j).read_index));
            ++i;++j;
        }
        else if(reads_on_SNP.at(now_snp_name).at(i).read_index < reads_on_SNP.at(next_snp_name).at(j).read_index){
            ++i;
        }
        else{
            ++j;
        }
    }
}

void Kernel::find_same_reads_in_accumlate(){
    forward_intersect_read.clear();
    unsigned int next_snp_name = positions->at(snp_iterator+1);
    size_t i =0;

    while(i < reads_on_SNP[next_snp_name].size()) {
        if(std::find(current_partition[0][0].begin(),current_partition[0][0].end(),reads_on_SNP[next_snp_name][i].read_index) != current_partition[0][0].end()){
            forward_intersect_read.push_back(reads_on_SNP[next_snp_name][i].read_index);
        }
        else if(std::find(current_partition[0][1].begin(),current_partition[0][1].end(),reads_on_SNP[next_snp_name][i].read_index) != current_partition[0][1].end()){
            forward_intersect_read.push_back(reads_on_SNP[next_snp_name][i].read_index);         
        }
        ++i;
    }
}

inline void Kernel::make_encode(const int &pos,const uint64_t &enco){
    encode_on_cur_partition.insert(pair<int,uint64_t>(pos,enco));
    return; 
}

inline bool Kernel::find_read_in_parition(int k , int iterator){
    if(k == best_answer[iterator]-1){
        return true;
    }
    /*if(best_answer[iterator] <= reads_on_SNP[positions->at(iterator)].size() ){
        return false;
    }*/
    /*if(k == partition_move_two.at(positions->at(iterator)).at(best_answer.at(iterator)-1-reads_on_SNP[positions->at(iterator)].size()).first
                    || k == partition_move_two.at(positions->at(iterator)).at(best_answer.at(iterator)-1-reads_on_SNP[positions->at(iterator)].size()).second)
        return true;*/
    return false;
}
inline bool Kernel::inverse_value(int pos){

    if(pos <= 0)
        return false ;

    if(trace_back_partition[pos-1][best_answer[pos]] == 1)
        return true;

    return false;
}
void Kernel::trace(){
    // accomplish block
    blockset.add_block(new Block(0,-1) , positions->at(snp_size-1) , block_size);
    // start at the end
    while(1){

        best_answer.push_back(best_pos);

        if(trace_back.size() == 0){
            break;
        }
        // pop the trace back vector
        auto record = trace_back[trace_back.size()-1];
        best_pos = record[best_pos];
        trace_back.pop_back();
    }

    std::reverse(best_answer.begin() , best_answer.end());

}

void Kernel::using_best_answer_separate_two_group() {
    
    for(unsigned int snp_id = 0;snp_id <snp_size -1;++snp_id){
        if(inverse_value(snp_id)){
            inverse_state = !inverse_state;
        }
        for(size_t read_id=0;read_id<reads_on_SNP.at(positions->at(snp_id)).size();++read_id){
            //cout<<read_id<<endl;
            if(__check_group(snp_id,read_id,1)){
                std::unordered_map <int, pair<int,int>>::iterator it = partition_read.find(reads_on_SNP.at(positions->at(snp_id)).at(read_id).read_index); 
                if (it != partition_read.end())
                    it->second.first ++;
                else
                partition_read.insert(pair<int, pair<int,int>>(reads_on_SNP.at(positions->at(snp_id)).at(read_id).read_index,pair<int,int>(1,0)));
            }
            else{
                std::unordered_map <int, pair<int,int>>::iterator it = partition_read.find(reads_on_SNP.at(positions->at(snp_id)).at(read_id).read_index); 
                if (it != partition_read.end())
                    it->second.second ++;
                else
                    partition_read.insert(pair<int, pair<int,int>>(reads_on_SNP.at(positions->at(snp_id)).at(read_id).read_index,pair<int,int>(0,1)));
            }
        }
    }

}
void Kernel::re_partition(){
    // parition again by mini movement
    
    for(auto &read:partition_read){
        
        int readID = read.first;
        int group_1 = read.second.first;
        int group_2 = read.second.second;
        //cout<<readID<<" G1: "<<group_1<<" G2: "<<group_2<<endl;
        if(group_1 >= group_2){
            readset.insert(pair<int, int>(readID,0));
        }
        else{
            readset.insert(pair<int, int>(readID,1));
        }
    }
    
}
void Kernel::determine_value_depend_on_group(int &value_a , int &value_b ,int snp_id){
    int count_A_0=0;
    int count_A_1=0;
    int count_B_0=0;
    int count_B_1=0;

    //#pragma omp for reduction( +:count_A_0,count_A_1,count_B_0,count_B_1)
    for(size_t k=0;k<reads_on_SNP.at(positions->at(snp_id)).size();++k){
        if((readset.at(reads_on_SNP.at(positions->at(snp_id)).at(k).read_index)) == 0){
            if(reads_on_SNP.at(positions->at(snp_id)).at(k).genotype == '0'){
                ++count_A_0;
            }
            else{
                ++count_A_1;
            }
        }
        else{
            if(reads_on_SNP.at(positions->at(snp_id)).at(k).genotype == '0'){
                ++count_B_0;
            }
            else{
                ++count_B_1;
            }
        }
        
    }

    if(count_A_0 + count_B_1 < count_B_0 + count_A_1){
        value_a = 1;
        value_b = 0;
    }
    else{
        value_b = 1;
        value_a = 0;
    }

    
}
void Kernel::get_read_info(unique_ptr<vector<const Entry *>>&input_column, int snp_iterator){
    vector<ReadInfo> readInfo_vec;
    // SNP position is a map index. map.at(index) is a read vector 
    reads_on_SNP.insert(pair<unsigned int , vector<ReadInfo>>(positions->at(snp_iterator) , readInfo_vec));
    for(size_t read_iter =0;read_iter<input_column->size();++read_iter){
        ReadInfo readInfo_tmp;
        readInfo_tmp.cost = input_column->at(read_iter)->get_phred_score();
        readInfo_tmp.read_index = input_column->at(read_iter)->get_read_id();
        readInfo_tmp.genotype = input_column->at(read_iter)->get_allele_type()+'0';
        // for now , this code just can handle phasing problem which only contain ALT & REF
        if((input_column->at(read_iter)->get_allele_type() == Entry::REF_ALLELE || input_column->at(read_iter)->get_allele_type() == Entry::ALT_ALLELE ))
            reads_on_SNP.at(positions->at(snp_iterator)).push_back(readInfo_tmp);

    }/*
    fstream file;
    file.open("gary_test/mutiple.txt", ios::app);
	file<<"input : \n";
    
    file<<"snp : "<<(positions->at(snp_iterator))<<endl;
    for(int i =0 ; i<reads_on_SNP.at(positions->at(snp_iterator)).size() ; ++i){
        file<<reads_on_SNP.at(positions->at(snp_iterator)).at(i).read_index<<" ";
    }
    file<<"\n";
    for(int i =0 ; i<reads_on_SNP.at(positions->at(snp_iterator)).size() ; ++i){
        file<<reads_on_SNP.at(positions->at(snp_iterator)).at(i).genotype<<" ";
    }
    file<<"\n";
    file.close();*/
}
void Kernel::pass_data(){
    previous_partition.clear();
    previous_partition = std::move(current_partition); 
    current_partition.clear();
    encode_on_pre_partition.clear();
    encode_on_pre_partition = std::move(encode_on_cur_partition); 
    encode_on_cur_partition.clear();
    previous_snp_movement.clear();
    previous_snp_movement = std::move(current_snp_movement); 
    current_snp_movement.clear();
    backward_intersect_read.clear();
    backward_intersect_read = std::move(forward_intersect_read); 
    forward_intersect_read.clear();
    snp_iterator++;
}

inline bool Kernel::__check_group(int snp_id , int read_id , int type){
    //cout<<"now check "<<reads_on_SNP.at(positions->at(snp_id)).at(read_id).read_index<<endl;
    
    if(type == 1){
        //cout<<(reads_on_SNP.at(positions->at(snp_id)).at(read_id).genotype == '0')<<" ";
        //cout<<find_read_in_parition( read_id , snp_id)<<" ";
        //cout<<inverse_state<<"\n";
        if(((reads_on_SNP.at(positions->at(snp_id)).at(read_id).genotype == '0') && !(find_read_in_parition( read_id , snp_id) ^ inverse_state))
                    ||  ((reads_on_SNP.at(positions->at(snp_id)).at(read_id).genotype == '1') && (find_read_in_parition( read_id , snp_id) ^ inverse_state) ))
            {
                return true;
            }
        return false ;
    }
    else{
        //cout<<(reads_on_SNP.at(positions->at(snp_id)).at(read_id).genotype == '1')<<" ";
        //cout<<find_read_in_parition( read_id , snp_id)<<" ";
        //cout<<inverse_state<<"\n";
        if(((reads_on_SNP.at(positions->at(snp_id)).at(read_id).genotype == '1') && !(find_read_in_parition( read_id , snp_id) ^ inverse_state))
                    ||  ((reads_on_SNP.at(positions->at(snp_id)).at(read_id).genotype == '0') && (find_read_in_parition( read_id , snp_id) ^ inverse_state) ))
            {
                return true;
            }
        return false ;
    }
    
}
void Kernel::insertionSort(vector<int> &data){
    //auto begin = std::chrono::high_resolution_clock::now();
    size_t i, j;
    int tmp;
    for(i = 1; i < data.size(); i++){
        tmp = data[i];
        for( j=i; j > 0 && tmp < data[j-1]; j-- )
            data[j] = data[j-1];        
        data[j] = tmp;
    }
    //auto end = std::chrono::high_resolution_clock::now();
    //auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    //cout << setiosflags(ios::left)<< setw(40) << "insertionSort: " <<setprecision(2) << elapsed.count()  <<"s"<<endl;
}