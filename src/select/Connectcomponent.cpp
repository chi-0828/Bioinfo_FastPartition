#include <iostream> 
#include <bits/stdc++.h>
#include <vector>
#include <assert.h>
#include <set>
#include "Connectcomponent.h"
int Connectcomponent::_find_node(map<int,Node> &nodes,int element){
    Node root = nodes[element];
    Node node = nodes[element];
    
    while(root.parent!=-1)
        root = nodes[root.parent];
            
    while(node.parent!=-1){
        node = nodes[node.parent];
        node.parent = nodes[root.node_num].parent;
    }
    return root.node_num;

}
void Connectcomponent::merge(map<int,Node> &node,int32_t first_node,int32_t second_node){
    assert( first_node!=second_node);

    int first_root = _find_node(node,first_node);
    int second_root = _find_node(node,second_node);
    if(first_root == second_root)
        return;

    if(node[first_node].node_num>node[second_node].node_num)
        node[first_node].parent = second_root;
       
    else{
        node[second_node].parent = first_root;
    }
}

void Connectcomponent::set_node(map<int,Node> &node,set<int> &positions){
    
    
    for(set<int>::iterator it = positions.begin();it!=positions.end();it++ ){
        Node insert;
        insert.node_num=*it;
        insert.parent=-1;
        node.insert(make_pair(*it,insert));
    }
}
int Connectcomponent::find(map<int,Node> &nodes,int id){
     int root = _find_node(nodes,id);
     return root;
}