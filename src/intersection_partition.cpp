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

Intersection_partition::Intersection_partition(ReadSet* read_set, const vector<unsigned int>* positions) :	
    read_set(read_set),
	input_column_iterator(*read_set, positions)
{
	read_set->reassignReadIds();

    std::vector<pair<unsigned int , string>> read_sourceID_Name;

	for (size_t i=0; i<read_set->size(); ++i) {
        read_sourceID_Name.push_back(pair<unsigned int , string>(read_set->get(i)->getSourceID() , read_set->get(i)->getName()));
	}

    cout<<"Get the final position of each read .... \n";
    for(auto read_id_name : read_sourceID_Name){
        auto Read = read_set->getByName(read_id_name.second,read_id_name.first);
        snp.read_last_pos.push_back(Read->lastPosition());
    }

    cout<<"Start algorithm\n";
	compute_table();
}
unique_ptr<vector<unsigned int> > Intersection_partition::extract_read_ids(const vector<const Entry *>& entries) {
	unique_ptr<vector<unsigned int> > read_ids(new vector<unsigned int>());
	for (size_t i=0; i<entries.size(); ++i) {
		read_ids->push_back(entries[i]->get_read_id());
	}
	return read_ids;
}

void Intersection_partition::compute_table() {
    // go to start point
	input_column_iterator.jump_to_column(0);
    snp.positions = input_column_iterator.get_positions();
    snp.snp_size = snp.positions->size();
    
    // calculate time cost
    double assenble_s=0;
    double partition_s=0;
    double build_s=0;
    double find_same_s=0;
    double DP_s=0;
    auto begin = std::chrono::high_resolution_clock::now();
    unique_ptr<vector<const Entry *> > current_input_column;
	unique_ptr<vector<const Entry *> > next_input_column;
	// get the next column ahead of time
	next_input_column = input_column_iterator.get_next();
	unique_ptr<vector<unsigned int> > next_read_ids = extract_read_ids(*next_input_column);
    // start our phasing algorithm. 
    // by GARY 2021-1-13
    while (snp.snp_iterator < snp.snp_size-1) {
        //cout<<"snp : "<<snp.positions->at(snp.snp_iterator)<<endl;
        // store info for debug (mutiple path)
        //map<int,vector<int>> temp_map;
        //snp.mutiple_choose.insert(pair<int,map<int,vector<int>>>(snp.snp_iterator,temp_map));
        ///////////////////////////////////////////////////////////////////////////////////////
        // current_input_column parameter will store each SNP and it mapped reads
        current_input_column = std::move(next_input_column);
        // get some information that we need.
        // i.e. read_id, read_genotype (0 or 1)
        if(snp.snp_iterator == 0){
            snp.get_read_info(current_input_column,snp.snp_iterator);
        }
        if (input_column_iterator.has_next()) {
            next_input_column = input_column_iterator.get_next();
		    next_read_ids = extract_read_ids(*next_input_column);
            /* Get next column read ID */
            snp.get_read_info(next_input_column,snp.snp_iterator+1);
        }
        //after first snp
        //do partition in current snp
        auto begin = std::chrono::high_resolution_clock::now();
        snp.partition();
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        partition_s += elapsed.count();

        //build encode
        begin = std::chrono::high_resolution_clock::now();
        if(snp.snp_iterator != 0){
            snp.build_table(PREVIOUS);
            end = std::chrono::high_resolution_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            build_s += elapsed.count();
            // combine two partition together, this parameter will store all appear read id.
            // In order to trace back to the best path in the future.
            begin = std::chrono::high_resolution_clock::now();
            snp.assemble();
            end = std::chrono::high_resolution_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            assenble_s += elapsed.count();

            // done and find best position
            // each partition have the cost of finally result 
            // minimum cost is the best result in the end
            if(!input_column_iterator.has_next()){
                // in order to find minimum cost, init a big value
                unsigned long long  min = snp.current_snp_movement.at(0);
                int best_pos_tmp = -1;
                for(auto const map:snp.current_snp_movement){
                    if(map.second >= 0)
                        if(map.second <= min){
                            min = std::move(map.second);
                            best_pos_tmp =std::move( map.first);
                        }
                }
                cout<<"best partition ID : "<<best_pos_tmp<<endl;
                cout<<"minimum  movement : "<<min<<endl;
                snp.best_pos = best_pos_tmp;
                break;
            }
        }
        //find same
        if(input_column_iterator.has_next()){
            begin = std::chrono::high_resolution_clock::now();
            if(snp.snp_iterator == 0){
                snp.find_same_reads();
            }
            else{
                snp.find_same_reads_in_accumlate();
            }
            end = std::chrono::high_resolution_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            find_same_s += elapsed.count();
        }
        //build encode for next round
        begin = std::chrono::high_resolution_clock::now();
        snp.build_table(NEXT);
        end = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        build_s += elapsed.count();

        snp.pass_data();

        //getchar();
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


std::pair<Read*,Read*> Intersection_partition::get_outputread() {
    auto tbegin = std::chrono::high_resolution_clock::now();
    cout<<"Get ouputread ..... \n";
    Read* Readset1 = new Read("ouputread_0", -1, -1,0);
    Read* Readset2 = new Read("ouputread_1", -1, -1,0);
	std::pair<Read*,Read*>ouputread(Readset1 , Readset2) ;

    // start our phasing algorithm-trace back. 
    // by GARY 2021-03-02  
    // find best answer for partition
    auto begin = std::chrono::high_resolution_clock::now();
    snp.trace();
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    cout << setiosflags(ios::left)<< setw(40) << "trace: " << fixed<<setprecision(2)  << elapsed.count() * 1e-9 <<"s"<<endl;
    // using answer to parition 
    // two group : readset1 and readset2
    begin = std::chrono::high_resolution_clock::now();
    snp.using_best_answer_separate_two_group();
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    cout << setiosflags(ios::left)<< setw(40) << "using_best_answer_separate_two_group: " << fixed<<setprecision(2)  << elapsed.count() * 1e-9 <<"s"<<endl;

    // depend on those error partition , redo it again 
    begin = std::chrono::high_resolution_clock::now();
    snp.re_partition();
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    cout << setiosflags(ios::left)<< setw(40) << "re_partition: " << fixed<<setprecision(2)  << elapsed.count() * 1e-9 <<"s"<<endl;


    // start write output string
    begin = std::chrono::high_resolution_clock::now();
    for(int i = 0;i <snp.snp_size -1;++i){

        int value_a = -1;
        int value_b = -1;
        // stroe value to superread
        snp.determine_value_depend_on_group(value_a ,value_b,i);
        if(value_a == -1 || value_b==-1 || (value_a==value_b)){
            cout<<value_a<<"!\n";
        }

		ouputread.first->addVariant((snp.positions->at(i)),value_a, 0);
		ouputread.second->addVariant((snp.positions->at(i)), value_b, 0);

    }

    // put into output readset

    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    cout << setiosflags(ios::left)<< setw(40) << "write: " << fixed<<setprecision(2)  << elapsed.count() * 1e-9 <<"s"<<endl;
    /*
    fstream file;
    file.open("gary_test/mutiple.txt", ios::app);
    file<<"answer : \n"<<output_read_set->at(0)->toString()<<endl;
    file.close();
    //print two group of reads

    file.open("gary_test/reads.txt", ios::out|ios::trunc);
    file<<"readset1 :"<<endl;
    for(int i=0;i<snp.readset1.size();i++){
        file<<snp.readset1.at(i)<<" ";
    }
    file<<"\n\nreadset2 :"<<endl;
    for(int i=0;i<snp.readset2.size();i++){
        file<<snp.readset2.at(i)<<" ";
    }
    file<<"\n\n";
    file.close();*/
    auto tend = std::chrono::high_resolution_clock::now();
    auto telapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(tend - tbegin);
    cout << setiosflags(ios::left)<< setw(40) << "trace back: " <<setprecision(2) << telapsed.count() * 1e-9 <<"s"<<endl;

    return ouputread;
    // end our algorithm-trace back.
}

