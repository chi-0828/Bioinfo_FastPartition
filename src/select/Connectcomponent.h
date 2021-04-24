#include <iostream> 
#include <bits/stdc++.h>
#include <vector>
#include <set>
#ifndef _CONNECTCOMPONENT_
#define _CONNECTCOMPONENT_
using namespace std;
typedef struct node{
int node_num;
int parent;

}Node;

class Connectcomponent{
    
    public: Connectcomponent(){
        ;
    }
    
    int _find_node(map<int,Node> &nodes,int element);
    void merge(map<int,Node> &nodes,int32_t first_node,int32_t second_node);
    void set_node(map<int,Node> &nodes,set<int> &positions);
    map<int,Node> nodes;
    int find(map<int,Node> &nodes,int id);

};
#endif