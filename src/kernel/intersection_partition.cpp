#include <stdexcept>
#include <cassert>
#include <limits>
#include <fstream>
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>
#include <list>  
#include <bitset>  
#include "intersection_partition.h"
#include "kernel.h"

Intersection_partition::Intersection_partition(vector<ReadVariant> * read_set , std::vector<unsigned int> read_last_pos ,std::vector<unsigned int> read_first_pos) :	
    read_set(read_set)
{
    auto begin = std::chrono::high_resolution_clock::now();
    kernel_table = new Kernel();
    kernel_table->read_last_pos = read_last_pos;
    kernel_table->read_first_pos = read_first_pos;
    std::set<int> all_snp_records;

    for(std::vector<ReadVariant>::iterator readIter = read_set->begin() ; readIter != read_set->end() ; readIter++ ){
        for(std::vector<Variant>::iterator variantIter = (*readIter).variantVec.begin(); variantIter != (*readIter).variantVec.end() ; variantIter++ ){
            //std::cout<<"intsert "<<variantIter->position<<endl;;
            all_snp_records.insert(variantIter->position);
            
            auto find_snp = kernel_table->reads_on_SNP.find(variantIter->position);
            if( find_snp != kernel_table->reads_on_SNP.end()){
                ReadInfo spot_remp;
                spot_remp.read_index = (*readIter).ID;
                spot_remp.genotype = variantIter->allele+'0';
                spot_remp.cost = variantIter->quality;
                kernel_table->reads_on_SNP.at(variantIter->position).push_back(spot_remp);
                //cout<<(*readIter).read<<" ";
            }
            else{
                vector<ReadInfo> read_temp;
                ReadInfo spot_remp;
                spot_remp.read_index = (*readIter).ID;
                spot_remp.genotype = variantIter->allele+'0';
                spot_remp.cost = variantIter->quality;
                read_temp.push_back(spot_remp);
                kernel_table->reads_on_SNP.insert(std::make_pair(variantIter->position , read_temp));
            }
            
        }
    }
    
    for(auto &snp: all_snp_records){
        kernel_table->positions->push_back(snp);
        //cout<<snp<<" ";
        /*auto a =  kernel_table->reads_on_SNP.find(snp);
        cout<< snp <<" ";
        for(auto i : (*a).second ){
            cout<<i.read_index <<" ";
        }
        cout <<"\n---------------\n";*/
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    cout << setiosflags(ios::left)<< setw(40) << "prepare data : " <<setprecision(2) << elapsed.count() * 1e-9 <<"s"<<endl;
    
    cout<<"Using "<<read_set->size()<<" reads and "<<all_snp_records.size()<<" snp ... \n";
    cout<<"Start algorithm ... \n";
	compute_table();
}
Intersection_partition::~Intersection_partition(){
    delete kernel_table;
    delete read_set;
}

void Intersection_partition::compute_table() {
    // get whole table size
    kernel_table->snp_size = kernel_table->positions->size();
    
    // calculate time cost
    double assenble_s=0;
    double partition_s=0;
    double build_s=0;
    double find_same_s=0;
    double DP_s=0;
    auto begin = std::chrono::high_resolution_clock::now();
    // start our phasing algorithm. 
    // by GARY 2021-1-13
    while (kernel_table->snp_iterator < kernel_table->snp_size-1) {

        //do partition in current snp
        auto begin = std::chrono::high_resolution_clock::now();
        kernel_table->partition();

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        partition_s += elapsed.count();
        
        begin = std::chrono::high_resolution_clock::now();
        // bulid block record 
        block_start.insert(make_pair(kernel_table->positions->at(kernel_table->snp_iterator) , kernel_table->block_num));
        if(kernel_table->snp_iterator == 0){
            kernel_table->blockset.add_block(new Block(kernel_table->positions->at(kernel_table->snp_iterator) , kernel_table->block_num++));
            kernel_table->block_size = 1;
        }
        else{
            //build encode
            kernel_table->build_table(PREVIOUS);
            end = std::chrono::high_resolution_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            build_s += elapsed.count();
            // combine two partition together, this parameter will store all appear read id.
            // In order to trace back to the best path in the future.
            begin = std::chrono::high_resolution_clock::now();
            kernel_table->assemble();
            end = std::chrono::high_resolution_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            assenble_s += elapsed.count();
      
            // done and find best position
            // each partition have the cost of finally result 
            // minimum cost is the best result in the end
            if(kernel_table->snp_iterator == kernel_table->snp_size-1-1){
                // in order to find minimum cost, init a big value
                unsigned long long  min = 1000000;
                int best_pos_tmp = -1;
                for(auto const map:kernel_table->current_snp_movement){
                    if(map.second <= min){
                        min = std::move(map.second);
                        best_pos_tmp =std::move( map.first);
                    }
                }
                cout<<"best partition ID : "<<best_pos_tmp<<endl;
                cout<<"minimum  movement : "<<min<<endl;
                kernel_table->best_pos = best_pos_tmp;
                break;
            }
        }
        //find same
        if(kernel_table->snp_iterator < kernel_table->snp_size-1-1){
            begin = std::chrono::high_resolution_clock::now();
            if(kernel_table->snp_iterator == 0){
                kernel_table->find_same_reads();
            }
            else{
                kernel_table->find_same_reads_in_accumlate();
            }
            end = std::chrono::high_resolution_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            find_same_s += elapsed.count();
        }
        //build encode for next round
        begin = std::chrono::high_resolution_clock::now();
        kernel_table->build_table(NEXT);
        end = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        build_s += elapsed.count();

        kernel_table->pass_data();
    }

    //print time cost
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    DP_s += elapsed.count();
    cout << setiosflags(ios::left)<< setw(40) << "DP loop: " << fixed<<setprecision(2)  << DP_s * 1e-9 <<"s"<<endl;
    cout << setiosflags(ios::left)<< setw(40) << "assemble: " << fixed<<setprecision(2)   << assenble_s * 1e-9 <<"s"<<endl;
    cout << setiosflags(ios::left)<< setw(40) << "build_encode: " << fixed<<setprecision(2)   << build_s * 1e-9 <<"s"<<endl;
    cout << setiosflags(ios::left)<< setw(40) << "find same: " << fixed<<setprecision(2)   << find_same_s * 1e-9 <<"s"<<endl;
    cout << setiosflags(ios::left)<< setw(40) << "parition: " << fixed<<setprecision(2)   << partition_s * 1e-9 <<"s"<<endl;
    
    // end our algorithm.
}


void Intersection_partition::get_outputread(std::string vcfFile) {
    auto tbegin = std::chrono::high_resolution_clock::now();
    cout<<"Get ouputread ..... \n";
    
    //map<pair<int, int>,int> result ;
    // start our phasing algorithm-trace back. 
    // by GARY 2021-03-02  
    // find best answer for partition
    kernel_table->trace();

    // using answer to parition 
    // two group : readset1 and readset2
    kernel_table->using_best_answer_separate_two_group();

    // depend on those error partition , redo it again 

    kernel_table->re_partition();

    // start write output string
    for(unsigned int i = 0;i <kernel_table->snp_size -1;++i){

        int value_a = -1;
        int value_b = -1;
        // stroe value to superread
        kernel_table->determine_value_depend_on_group(value_a ,value_b,i);
        if(value_a == -1 || value_b==-1 || (value_a==value_b)){
            cout<<value_a<<"!\n";
        }
        unsigned int pos = kernel_table->positions->at(i);

        auto pos_allele = make_pair(value_a,value_b);

        result.insert(make_pair(pos,pos_allele));
    }

    auto tend = std::chrono::high_resolution_clock::now();
    auto telapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(tend - tbegin);
    cout << setiosflags(ios::left)<< setw(40) << "trace back: " <<setprecision(2) << telapsed.count() * 1e-9 <<"s"<<endl;

    cout<<"blocks num : "<<kernel_table->blockset.get_size()<<endl;
    cout<<kernel_table->blockset.find_largest_block()<<endl;;
    // end our algorithm-trace back.

    auto wbegin = std::chrono::high_resolution_clock::now();
    writingResult(vcfFile);
    auto wend = std::chrono::high_resolution_clock::now();
    auto welapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(wend - wbegin);
    cout << setiosflags(ios::left)<< setw(40) << "write VCF: " <<setprecision(2) << welapsed.count() * 1e-9 <<"s"<<endl;
}

void Intersection_partition::writingResult(std::string vcfFile){

    std::ifstream file(vcfFile);
    std::ofstream resultVcf("result.vcf");

    if(!file.is_open()){
        std::cout<< "Fail to open file: " << vcfFile << "\n";
    }
    else if(!resultVcf.is_open()){
        std::cout<< "Fail to open vcf: result.vcf\n";
    }
    else{
        std::string input;
        while(! file.eof() ){
            std::getline(file, input);
            
            if( input.find("#")!= std::string::npos){
                // header
                if(input.find("##")== std::string::npos){
                    resultVcf <<  "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">\n";
                }
                resultVcf << input << "\n";
            }
            else if( input.find("#")== std::string::npos ){
                std::istringstream iss(input);
                std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

                if( fields.size() == 0 )
                    break;

                int pos = std::stoi( fields[1] );


                fields[8] = fields[8] + ":PS";
                std::map<int , pair<uint8_t, uint8_t>>::iterator Iter = result.find(pos-1);
                auto block_pos = block_start.find(pos);
                if( Iter != result.end() ){
                    fields[9] = fields[9] + ":" + std::to_string((*block_pos).second);
                    
                     // find GT flag
                    int colon_pos = 0;
                    int gt_pos = fields[8].find("GT");
                    for(int i =0 ; i< gt_pos ; i++){
                        if(fields[8][i]==':')
                            colon_pos++;
                    }
                    // find GT value start
                    int current_colon = 0;
                    int modifu_start = 0;
                    for(unsigned int i =0; i < fields[9].length() ; i++){
                        if( current_colon >= colon_pos )
                            break;
                        if(fields[9][i]==':')
                            current_colon++;  
                        modifu_start++;
                    }

                    fields[9][modifu_start] =  (Iter->second.first) +'0';
                    fields[9][modifu_start+1] = '|';
                    fields[9][modifu_start+2] =  (Iter->second.second) +'0';
                    
                    
                    //std::cout<< pos << "\t" << subNodeHP[tmpA] << "|" << subNodeHP[tmpB] << "\n";
                    for(std::vector<std::string>::iterator fieldIter = fields.begin(); fieldIter != fields.end(); ++fieldIter){
                        if( fieldIter != fields.begin() )
                            resultVcf<< "\t";
                        resultVcf<< (*fieldIter);
                    }
                    resultVcf<< "\n";
                }
                else{
                    resultVcf<< input << "\n";
                }
            }
        }
    }
}