#include <iostream> 
#include <bits/stdc++.h>
#include <vector>
#include <set>
#include <map>
#include "../phasing/Util.h"

using namespace std;

#ifndef _PRIORITYQUEUE_
#define _PRIORITYQUEUE_


typedef int item_type;
typedef vector<int> priority_type;
//typedef pair<priority_type_ptr,item_type> queue_entry_type;

typedef priority_type* priority_type_ptr;

class Priorityqueue{
	
	public:
		vector<pair<priority_type_ptr,item_type>> heap;
		unordered_map<item_type,int> positions;
		void c_push(priority_type_ptr score, int item);
		void _swap(int index1, int index2);
		bool _score_lower( int index1, int index2);
		void _sift_up( int index);
		void _sift_down( int index);
		void c_change_score(item_type &item, priority_type_ptr &c_new_score);
		pair<priority_type_ptr,item_type> c_pop();
		priority_type_ptr c_get_score_by_item(item_type item);
		int size(); 
		bool c_is_empty();
		bool is_empty();
		void __dealloc__();


};
#endif 