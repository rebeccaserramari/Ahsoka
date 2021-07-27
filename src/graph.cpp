#include "graph.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <vector>
#include <utility>
#include <set> 
#include <unordered_set> 

using namespace std;


DirectedNode::DirectedNode() :
id(0),
end(false)
{
}

DirectedNode::DirectedNode(int id, bool end) :
id(id),
end(end)
{
}

Node::Node():
node_id(0),
bubble_id(0),
chain_id(0),
allele(0),
node_seq(""),
childrenleft(),
childrenright(),
visited(false)
{
}

Node::Node(int id, string seq):
node_id(id),
bubble_id(0),
chain_id(0),
allele(0),
node_seq(seq),
childrenleft(),
childrenright(),
visited(false)
{
}

bool Node::operator==(const Node& other) const
{
	return (node_id == other.node_id && node_seq == other.node_seq);
}

bool Node::operator!=(const Node& other) const
{
	return !(*this == other);
}

bool Node::operator<(const Node& other) const
{
	return(node_id < other.node_id);
}

bool Node::operator>(const Node& other) const
{
	return(node_id > other.node_id);
}

bool DirectedNode::operator==(const DirectedNode& other) const
{
	return(id == other.id && end == other.end);
}

bool DirectedNode::operator!=(const DirectedNode& other) const
{
	return !(*this == other);
}


Bubble::Bubble():
	id(0),
	source(),
	sink(),
	innerNodes()
{
}

Bubble::Bubble(Node source, Node sink, std::vector<Node> innerNodes):
	id(0),
	source(source),
	sink(sink),
	innerNodes(innerNodes)
	{
	}

vector<Node> Bubble::getNodes(){
	vector<Node> nodes;	
	nodes.push_back(source);
	for (auto& node: innerNodes)
		nodes.push_back(node);
	nodes.push_back(sink);
	return(nodes);
}

vector<int> Bubble::getNodeIds() {
	vector<int> nodeIds;
	nodeIds.push_back(source.node_id);
	for (auto& node: innerNodes)
		nodeIds.push_back(node.node_id);
	nodeIds.push_back(sink.node_id);
	return(nodeIds);
}

Node Bubble::getNode(int node_id) {
	if (source.node_id == node_id)
		return(source);
	else if (sink.node_id == node_id)
		return(sink);
	else {
		for (auto& node: innerNodes){
			if (node.node_id == node_id)
				return(node);		
		}	
	}
}

Chain::Chain():
	id(0),
	bubbles()
	{
	}

void Chain::addBubble(Bubble b){
	bubbles.push_back(b);
}

vector<Node> Chain::listChain() {
	vector<Node> cnodes;	
	for (auto b: bubbles){
		for (auto n: b.getNodes())
			if (find(cnodes.begin(),cnodes.end(),n) == cnodes.end())
				cnodes.push_back(n);	
	}
	return(cnodes);
}

int Chain::size() {
	return(bubbles.size());
}

Graph::Graph():
	nodes(),
	edges(),
	offsets(),
	chains()
{
}

int raw_id(string s) {
//remove all non-numeric characters from s
s.erase(remove_if(s.begin(), s.end(), [](char c) { return !isdigit(c); }), s.end());
return stoi(s);
}

Chain Graph::getChain(int chain_id) {
	for (auto& chain: chains) {
		if (chain.id == chain_id)
			return(chain);	
	}
}

void Graph::printChain(int chain_id) {
	cout << "chain id: " << chain_id << endl;
	for (auto& bubble: getChain(chain_id).bubbles) {
		cout << "bubble id: " << bubble.id << endl;
		cout << "node id: " ;
		for (auto& node: bubble.getNodes())	{
			cout << node.node_id << ",";
		}
		cout << endl;
	}	
}



Graph Graph::ReadGraph(string filename) {
	ifstream file(filename);
	Graph graph;
	string line;
	while (getline(file, line)){
		
		if (!file.good()) break;
		if (line.size() == 0) continue;
		if (line[0] != 'S' && line[0] != 'L') continue;
		if (line[0] == 'S'){
			stringstream sstr {line};
			string idstr;
			string s;
			string seq;
			sstr >> s;
			assert(s == "S");
			sstr >> idstr;
			int id = raw_id(idstr);
			sstr >> seq;
			assert(seq.size() >= 1);

			//set id to node object and fill with sequence; later with bubble/chain/allele info
			graph.nodes[id] = Node(id, seq);
		}
		if (line[0] == 'L'){
			stringstream sstr {line};
			string l;
			string start;
			string start_order;
			string end;
			string end_order;
			int offset;
			sstr >> l;
			assert(l == "L");
			sstr >> start;
			int start_id = raw_id(start);
			sstr >> start_order;
			sstr >> end;
			int end_id = raw_id(end);
			sstr >> end_order;
			assert(start_order == "+" || start_order == "-");
			assert(end_order == "+" || end_order == "-");
			sstr >> offset;
			char dummyc;
			sstr >> dummyc;
			assert(dummyc == 'M' || (dummyc == 'S' && offset == 0));
			
			assert(offset >= 0);
			
			DirectedNode from(start_id, start_order == "+");
			DirectedNode to(end_id, end_order == "+");
			graph.edges[from].push_back(to);
			if (start_order == "+")
				graph.nodes[start_id].childrenleft.push_back(make_pair(graph.nodes[end_id].node_id, end_order == "+"));
			if (start_order == "-")
			//assert that end_id should == node_id
				graph.nodes[start_id].childrenright.push_back(make_pair(graph.nodes[end_id].node_id,end_order == "+"));
			graph.offsets[make_pair(from, to)] = offset;			
		}
	}
	return graph;
}

pair<DirectedNode,DirectedNode> Graph::getEdge(Node first, Node second) {
	for (auto& val: {true, false}) {	
		DirectedNode from(first.node_id, val) ;
		for (auto& to_node: edges[from]) {
			if (to_node.id == second.node_id){
				return(make_pair(from, to_node));		
			}	
		}
	}
	//return(0);
}

/*
def find_bubblechains(graphfile, mygraph):
	#TODO: k value is not correct here
	nodes = read_gfa(graphfile, 61)
	graph = Graph(graphfile)
	graph.nodes = nodes
	print("no. of nodes: ", len(nodes))
	print("no. of nodes in mygraph: ", len(mygraph.nodes))
	
	graph.find_chains()
	print("no. of bchains: ", len(graph.b_chains))
	nodecount = 0
	nodes_in_chains = []
	chains_to_end = {}
	endnodes = set()
	chainid = 1
	for chain in graph.b_chains:
		chain.set_id(chainid)
			
		#bubble chains need some identifier too? In the fasta they are labeled with chain0 to chain27		
		print("number of bubbles in chain: ", chainid, len(chain.sorted))
		nodecount += len(chain.list_chain())
		nodes_in_chains.extend(chain.list_chain())
		node_cn = {}
		chains_to_end[chain] = []
		bubbleid = 0
#		for bubble in chain.bubbles:
		for bubble in chain.sorted:	
			bubble.set_id(bubbleid)
			
			
			for n in [bubble.source]+[bubble.sink]+bubble.inside:
				
				mynode = mygraph.get_node(str(n.id))
				
				mynode.set_chainid(chainid)
				mynode.set_bubbleid(bubbleid)
				
			if (bubble.source.id not in node_cn.keys()):
				node_cn[bubble.source.id] = 1
			else:
				node_cn[bubble.source.id] += 1
			if (bubble.sink.id not in node_cn.keys()):
				node_cn[bubble.sink.id] = 1			
			else:			
				node_cn[bubble.sink.id] += 1
			bubbleid += 1
		chainid += 1
#		print("number of bubbles in chain: ", len(chain.bubbles))
#		print("number of node cn keys: ", len(node_cn.keys()))
		endcounter = 0	
		#those bubble sources and sinks that appear only once are the end of a bubble chain
		for end in node_cn.keys():
			if node_cn[end] ==1:
				chains_to_end[chain].append(end)
				buf = 6-len(str(end))
				endnodes.add("utg"+buf*"0"+str(end)+"l")
				endcounter += 1
			 
#		print("endcounter: ", endcounter)
	print("first 5 nodes in chains: ", nodes_in_chains[:5])			
	print("first 5 nodes in mygraph: ", [node.id for node in mygraph.nodes[:5]])
	print("first 5 nodes in mygraph: ", [node.raw_id for node in mygraph.nodes[:5]])
	print("number of keys in chains_to_end (should equal no. of bchains): ", len(chains_to_end.keys()))	
	print("number of nodes present in chains: ", nodecount)
	
	#returns all nodes that are present in bubble chains, dictionary mapping chain object to it's end nodes
	return(nodes_in_chains, chains_to_end, graph) 
*/

void Graph::findBubbles(){
	//find bchains
	for (auto node: nodes) {
	//	cout << "node.first: " << node.first << ", node.second.id: " << node.second.node_id << endl;
		Chain bchain;
		if (!node.second.visited) {
			vector<int> dirs = {0,1};
			for (auto i: dirs) {
				findBubble(node.second, i==1, &bchain);			
			}		
			if (bchain.size() != 0) {
				
				chains.push_back(bchain);
				for (auto n: bchain.listChain())
					n.visited = true;			
			}
		}	
	}
	int chain_id = 0;
	int bubble_id = 0;
	for (auto& chain: chains) {
		chain.id = chain_id;
		bubble_id = 0;
		for (auto& bubble: chain.bubbles) {
			bubble.id = bubble_id;
			for (auto& node: bubble.getNodes()) {
				nodes[node.node_id].chain_id = chain_id;			
				nodes[node.node_id].bubble_id = bubble_id;			
			}		
			bubble_id += 1;
		}
		chain_id += 1;	
	}
	
	//TODO add attr 'end' to bchains to denote bubble sources and sinks that mark the two ends of a chain (should be two each time?)
}

void Graph::findBubble(Node node, bool direction, Chain* bchain) {
	unordered_set<DirectedNode> seen;
	unordered_set<int> visited;
	vector<Node> nodesInside;
	seen.insert(DirectedNode(node.node_id, direction));
	set<pair<Node, bool>> S;
	S.insert(make_pair(node, direction));
	while (S.size() > 0){

		auto v = *S.begin();
		S.erase(S.begin());
		v = make_pair(v.first, v.second);
		visited.insert(v.first.node_id);
		nodes[v.first.node_id].visited = true;
		nodesInside.push_back(v.first);

		//node visited, so remove from the list of seen nodes
		DirectedNode nodepos_v = DirectedNode(v.first.node_id, v.second);
		seen.erase(nodepos_v);
		vector<pair<int, int>> children;
		//node children
		if (v.second == false)
			children = v.first.childrenleft;
		else
			children = v.first.childrenright;
		//in case of a tip
		if (children.size() == 0)
			break;
	//	cout << "v.first.id: " << v.first.node_id << ", v.second: " << v.second << endl;
		for (auto u: children) {
				vector<pair<int, int>> u_parents;
				int u_child_direction;
				if (u.second == 0){
					u_child_direction = 1;
					u_parents = nodes[u.first].childrenleft;
				}
				else{
					u_child_direction = 0;
					u_parents = nodes[u.first].childrenright;
				}
		//		cout << "u: " << u.first << endl;
		//		cout << "u child direction: " << u_child_direction << endl;
		//		cout << "u's parents: ";
		//		for (auto& p: u_parents) {
		//			cout << p.first << ",";				
		//		}
		//		cout << endl;
				if (u.first == node.node_id) {
					//loop found
					
					set<pair<Node, bool>> newS;
					S = newS;
					break;
				}
				
				if (u.second == 0)
					seen.insert(DirectedNode(u.first, true));
				else
					seen.insert(DirectedNode(u.first, false));
				//all u_parents are visited-> push into S
				bool all_visited = true;
				for (auto p: u_parents) {
		//			cout << "p: " << p.first << endl;
						//if one of the nodes is not visited
						// i is a pair(node_id, Node)
						if (visited.find(p.first) == visited.end())		
							all_visited = false;			
					
				}
				//exception for when there is one parent with two children of which one is a tip?
				
				
				bool single_parent_has_tip = false;
				vector<pair<int,int>> v_c;

				vector<pair<int, int>> v_parents;
				
				if (v.second == 0){
					v_parents = v.first.childrenleft;
				}				
				else{
					v_parents = v.first.childrenright;				
				}
				if (v_parents.size() == 1) {
				//	cout << "v: " << v.first.node_id << endl;
				//	for (auto& vp: v_parents)		{				
				//		cout << "parents, vp.first "<< vp.first << endl;
				//} 
					pair<int,int> parent = *v_parents.begin();
					if (v.second == false)
						v_c = nodes[parent.first].childrenleft;
					else
						v_c = nodes[parent.first].childrenright;
					//if v_c has 2 children, one should be v, the other a tip
					for (auto& i: v_c){
						if (nodes[i.first].childrenleft.size() == 0 || nodes[i.first].childrenright.size() == 0 )
							single_parent_has_tip = true;					
					}	
				}
			//TODO the line below strangely leads to erroneous bubble chains: the chains get partially listed double/one chain ends too early
			//	if (single_parent_has_tip) all_visited = false;
				if (all_visited)
					S.insert(make_pair(nodes[u.first], u_child_direction==1));

		}
		
		if ((S.size() == 1) && (seen.size() == 1)){
				auto t = *S.begin();
				S.erase(S.begin());
				nodesInside.push_back(t.first);
				
				if (nodesInside.size() == 2)
					break;

				vector<Node>::iterator node_pos = find(nodesInside.begin(), nodesInside.end(), node);
				if (node_pos != nodesInside.end())
    				nodesInside.erase(node_pos);
				vector<Node>::iterator t0_pos = find(nodesInside.begin(), nodesInside.end(), t.first);
				if (t0_pos != nodesInside.end())
    				nodesInside.erase(t0_pos);

				//create Bubble object with source, sink and inside)
				Bubble bubble(node, t.first, nodesInside);
				//maybe check if bubble already contained in chain	
							
				bchain->addBubble(bubble);

				findBubble(t.first, t.second, bchain);
		}
	}

}
Node Graph::getNode(int node_id) {
	bool node_found = false;	
	for (auto& it: nodes) {
		if (node_id == it.first) {
			node_found = true;
			return it.second;	
		}
	}
	if (!node_found) 
		std::cerr << "Node not in graph!" << std::endl;
}

void Graph::addNode(Node node) {
	int id = node.node_id;
	nodes[id] = node;
}