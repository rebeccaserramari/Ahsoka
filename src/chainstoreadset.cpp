#include <iostream>
#include <sstream>

#include "argumentparser.hpp"
#include "timer.hpp"

#include "graph.hpp"
#include <vector>
#include <unordered_map>
#include <algorithm>

using namespace std;

//TODO add alleles/ dict(allele, allelepath) to bubbles



unordered_map<int,vector<int>> findPathsSimple(Bubble bubble) {
	
	unordered_map<int,vector<int>> idPaths;
	
	Node node;
	for(int i = 0; i < bubble.innerNodes.size(); ++i) {
		node = bubble.innerNodes[i];
		idPaths[i].push_back(bubble.source.node_id);
		idPaths[i].push_back(node.node_id);
		idPaths[i].push_back(bubble.sink.node_id);	
	}  
	
	return(idPaths);
}


template <typename T>
bool contains(vector<T> vec, const T & elem)
{
    bool result = false;
    if( find(vec.begin(), vec.end(), elem) != vec.end() )
    {
        result = true;
    }
    return result;
}

void addSequence(pair<Node,int> sinkPair, Bubble bubble, vector<int>* idSeq, vector<vector<int>>* idPaths) {
	int direction = sinkPair.second;
	Node node = sinkPair.first; 
		//if node_id is not yet in vector

		if (find(idSeq->begin(), idSeq->end(), node.node_id) == idSeq->end()){
			idSeq->push_back(node.node_id);			
		}
		bool withinBubble = true;
		vector<pair<int, int>> nextnodes;
		if (direction == 0) {
			nextnodes = node.childrenright;
		}
		else {
			nextnodes = node.childrenleft;	
		}
		for (auto& child: nextnodes) {
			
			//if child node is not within the bubble nodes	
			if (!contains(bubble.getNodeIds(),child.first)) {
				withinBubble = false;
			}
				
		}
		if (nextnodes.size() > 0 && withinBubble) {
			for (auto& child: nextnodes) {	
    			int index = find(idSeq->begin(), idSeq->end(), node.node_id) - idSeq->begin();
    			vector<int> vec(index + 1);
				copy(idSeq->begin(), idSeq->begin() + index + 1, vec.begin());	
				pair<Node, int> childPair = make_pair(bubble.getNode(child.first), child.second);
				addSequence(childPair, bubble, &vec, idPaths); 			
			}		
		}
		else {
			idPaths->push_back(*idSeq);		
		}
				
	
}

vector<vector<int>> findPathsComplex(Bubble bubble) {
	vector<vector<int>> idPaths;
	//TODO is this correct? start from sink and only look at 0 (childrenright)?
	int direction = -1;
	bool right_withinbubble = true;
	bool left_withinbubble = true;
	for (auto& rchild: bubble.sink.childrenright) {
		if (!contains(bubble.getNodeIds(),rchild.first)) {
				right_withinbubble = false;
			}	
	}
	for (auto& lchild: bubble.sink.childrenleft) {
		if (!contains(bubble.getNodeIds(),lchild.first)) {
				left_withinbubble = false;
			}	
	}
	//if the right children are not within the bubble, go to the left
	if (!right_withinbubble) {
		direction = 1;	
	}
	//otherwise go to the right
	else direction = 0;

	pair<Node,int> sinkPair = make_pair(bubble.sink, direction);
	vector<int> idSeq;
	addSequence(sinkPair, bubble, &idSeq, &idPaths);
/*	if (idPaths[0].size() == 1) {
			pair<Node,int> sinkPair = make_pair(bubble.sink, 1);
			vector<int> idSeq;
			addSequence(sinkPair, bubble, &idSeq, &idPaths);
	}*/
	return(idPaths);
}


unordered_map<int, vector<vector<int>>> ChainsToReadset(Graph graph ) {
	int count = 0;
	unordered_map<int,vector<int>> idPaths;
	unordered_map<int, vector<vector<int>>> pathToAlleles;
	for (auto& chain: graph.chains) {
		count +=1 ;
		vector<int> seq;
		cout << "chain id: "<< chain.id << endl;
		for (auto& bubble: chain.bubbles ){
			
			vector<Bubble> bubblepaths;
			if (bubble.innerNodes.size() == 2) {
			idPaths = findPathsSimple(bubble);
			int allele = 0;
			vector<int> path;
				for (int i = 0; i < idPaths.size(); ++i) {
					path = idPaths[i];
					pathToAlleles[chain.id].push_back(path);	
					for (auto& nodeid: path) 
						graph.nodes[nodeid].allele = allele;				
					allele += 1;		
				}				
			}
			else {
			vector<vector<int>> idPathsComplex = findPathsComplex(bubble);	
			int allele = 0;
			vector<int> complexpath;
				for (int i = 0; i < idPathsComplex.size(); ++i) {
					complexpath = idPathsComplex[i];
					pathToAlleles[chain.id].push_back(complexpath);	
					for (auto& nodeid: complexpath) {
						graph.nodes[nodeid].allele = allele;	
							
					}		
					allele += 1;		
				}		
			}		
		}	
	}
	return(pathToAlleles);
}

unordered_map<int, unordered_map<int,vector<vector<int>>>> ChainsToReadsetDetailed(Graph graph ) {
	int count = 0;
	unordered_map<int,vector<int>> idPaths;
	vector<vector<int>> idPathsComplex;
	unordered_map<int, vector<vector<int>>> bubblepathToAlleles;
	unordered_map<int, unordered_map<int,vector<vector<int>>>> pathToAlleles;
	for (auto& chain: graph.chains) {
		count +=1 ;
		vector<int> seq;
		unordered_map<int, vector<vector<int>>> bubblepathToAlleles;
		for (auto& bubble: chain.bubbles ){
			if (bubble.innerNodes.size() == 2) {
			idPaths = findPathsSimple(bubble);
			int allele = 0;
			vector<int> path;
				for (int i = 0; i < idPaths.size(); ++i) {
					path = idPaths[i];
					bubblepathToAlleles[bubble.id].push_back(path);	
					for (auto& nodeid: path) 
						graph.nodes[nodeid].allele = allele;				
					allele += 1;		
				}
			pathToAlleles[chain.id] = bubblepathToAlleles;					
			}
			else {
			idPathsComplex = findPathsComplex(bubble);	
			int allele = 0;
			vector<int> complexpath;
				for (int i = 0; i < idPathsComplex.size(); ++i) {
					complexpath = idPathsComplex[i];
					bubblepathToAlleles[bubble.id].push_back(complexpath);	
					for (auto& nodeid: complexpath) {
						graph.nodes[nodeid].allele = allele;	
							
					}		
					allele += 1;		
				}
			pathToAlleles[chain.id] = bubblepathToAlleles;		
			}		
		}	
	}
	return(pathToAlleles);
}



