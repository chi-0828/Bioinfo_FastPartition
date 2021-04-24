#include <iostream> 
#include <bits/stdc++.h>
#include <vector>
#include <set>
#include <map>
#include <assert.h> 
#include "Priorityqueue.h"

using namespace std;

bool _vector_score_lower(priority_type_ptr &first, priority_type_ptr &second){
	for(long unsigned int i=0; i<min(first[0].size(), second[0].size()); i++){
		
		if (first[0][i] < second[0][i])
			return true;
		if (first[0][i] > second[0][i])
			return false;
	if (first[0].size() < second[0].size())
		return true;
	else
		return false;
    }
}
/*priority_type_ptr _pyscore_to_vector(priority_type_ptr score){
	
	priority_type_ptr result = new priority_type();
	if(isinstance(score,int))
		result[0].push_back(score);
	else{
            
			for(vector<int>::iterator it = score[0].begin;it!=score[0].end;it++)
				if(!isinstance(i,int))
					;
				result[0].push_back(i);
    }    
	return result;
    

}*/
int _parent(int index){
	//returns the index of the parent node in the heap depending on the given index
	if(index==0){
		return -1;
	}
	return (index - 1)/2;

}

int _left_child(int index){
	//return the index of the left child in the heap depending on the given index'''
	return (2 * index) + 1;
}


int _right_child(int index){
	//return the index of the right child in the heap depending on the given index'''
	return (2 * index) + 2;
}

//TODO: in general, see https://www.python.org/dev/peps/pep-0008

void Priorityqueue::__dealloc__(){
		for (long unsigned int i=0;i<this->heap.size();i++)
			delete this->heap[i].first;
			
}
void Priorityqueue::c_push(priority_type_ptr score, int item){
		//Pointer ownership is transferred to priorityqueue.
		int newindex = this->heap.size();
		pair<priority_type_ptr,item_type> entry;
		//cerr<<item;
		entry.first = score;
		entry.second = item ;
		this->heap.push_back(entry);
		this->positions[item] = newindex;
		this->_sift_up(newindex);
}
void Priorityqueue::_swap(int index1, int index2){
		//swaps the position of the nodes in the priority queue from first index(index1)and second index(index2)

		pair<priority_type_ptr,item_type> entry1 = this->heap[index1];
		int pos1 = this->positions[entry1.second];

		pair<priority_type_ptr,item_type> entry2 = this->heap[index2];
		int pos2 = this->positions[entry2.second];
		this->positions[entry1.second] = pos2;
		this->positions[entry2.second] = pos1;
		this->heap[index1] = entry2;
		this->heap[index2] = entry1;
}
bool Priorityqueue::_score_lower( int index1, int index2){
	    pair<priority_type_ptr,item_type> entry1 = this->heap[index1];
		pair<priority_type_ptr,item_type> entry2 = this->heap[index2];
		return  _vector_score_lower(entry1.first, entry2.first);
}
void Priorityqueue::_sift_up( int index){
		/*recursive method to check if score of item at given index is higher than the score of the parent
		and swaps the nodes till priorityqueue property is restored*/
		//print(str(this->), '_sift_up', index)
		int parentindex = _parent(index);
		//cerr<<parentindex;
		assert(parentindex !=index);
		if (parentindex >= 0){
			if(this->_score_lower(parentindex, index)){
				this->_swap(parentindex, index);
				this->_sift_up(parentindex);
			}
		}
}

void Priorityqueue::_sift_down(int index){
		/*Check if score of item at given index is lower than the score of its children,
		 therefore need to swap position with its children*/
		int rchildindex = _right_child(index);
		int lchildindex = _left_child(index);
		assert (rchildindex != index);
		assert (lchildindex != index);
		//if both children are in the heap
		//only need to know if right child exists then automatically also the left child exists
		if (rchildindex < this->heap.size()){
			if (this->_score_lower(lchildindex, rchildindex)){
				if (this->_score_lower(index, rchildindex)){
					this->_swap(rchildindex, index);
					this->_sift_down(rchildindex);
				}
			}
			else{
				if (this->_score_lower(index, lchildindex)){
					this->_swap(lchildindex, index);
					this->_sift_down(lchildindex);
				}
			}
		}
		else if( lchildindex < this->heap.size()){
			if (this->_score_lower(index, lchildindex)){
				this->_swap(lchildindex, index);
				this->_sift_down(lchildindex);
			}
		}
}
/*pair<priority_type_ptr,item_type> Priorityqueue::pop(){
		//Removes the item with largest score and returns it as tupel of (score, item).
		pair<priority_type_ptr,item_type> entry = this->c_pop();
		priority_type_ptr score = entry.first;
		item_type item = entry.second;
		if (score[0].size() == 1){
			result = (score[0][0], item)
		}
			
		else:
			result = tuple(score[0]), item
		delete score;
		return result
}*/
pair<priority_type_ptr,item_type> Priorityqueue::c_pop(){
		
		//Removes the item with largest score and returns it as tupel of (score, item).'''
		if(this->heap.size() == 0)
			cerr<<"PriorityQueue empty.";
		pair<priority_type_ptr,item_type> last_entry = this->heap[this->heap.size()-1];
		pair<priority_type_ptr,item_type> first_entry = this->heap[0];
		//looks if heap is only one element then no need restore heap property;
		if(this->heap.size() == 1){
			this->positions.erase(first_entry.second);
			this->heap.pop_back();
		}
		else{
			this->heap[0] = last_entry;
			this->heap.pop_back();

			this->positions[last_entry.second]= 0;
			this->positions.erase(first_entry.second);

			this->_sift_down(0);
		}
		return first_entry;
}

void Priorityqueue::c_change_score(item_type &item, priority_type_ptr &c_new_score){
	
		int position =  this->positions[item];
		priority_type_ptr c_old_score = this->heap[position].first;
		this->heap[position].first = c_new_score;

		//Differentiate between increasing and decreasing score
		if (_vector_score_lower(c_old_score, c_new_score))
			this->_sift_up(position);
		else
			this->_sift_down(position);

		delete c_old_score;
}

priority_type_ptr  Priorityqueue::c_get_score_by_item( item_type item){
		/*Returns score of the given item or NULL if item is not in the heap. 
		Pointer ownership stays with priorityqueue.*/
		unordered_map<item_type,int>::iterator it = this->positions.find(item);
		if (it == this->positions.end())
			return NULL;
		pair<priority_type_ptr,item_type> entry = this->heap[this->positions[item]];
		return entry.first;
		//priority_type_ptr score = entry.first;
		//return score;
		

}
int Priorityqueue::size(){
		return this->heap.size();
}
bool Priorityqueue::is_empty(){
	//Return if actual Priority Queue is Empty'''
	if(this->heap.size()==0)
		return true;
	else
		return false;
}

bool Priorityqueue::c_is_empty(){
	if(this->heap.size()==0)
		return true;
	else
		return false;
		
}
