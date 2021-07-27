#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <string>
#include <vector>
#include <algorithm>
#include <istream>
#include <ostream>
#include <unordered_map>


/** 
* Represents a Graph object
**/




class DirectedNode{
public:
	DirectedNode();
	DirectedNode(int id, bool end);
	int id;
	bool end;
	DirectedNode Reverse() const;
	bool operator==(const DirectedNode& other) const;
	bool operator!=(const DirectedNode& other) const;
};

template<>
struct std::hash<DirectedNode> {
		std::size_t operator()(const DirectedNode& x) const{
			return std::hash<int>()(x.id) ^ std::hash<bool>()(x.end);
		}
	};
	
template <> 
struct std::hash<std::pair<DirectedNode, DirectedNode>>
	{
		std::size_t operator()(const std::pair<DirectedNode, DirectedNode>& x) const
		{
			// simple hashing with hash<DirectedNode>()(x.first) ^ hash<DirectedNode>()(x.second) collides each edge formed like (x -> x+1)
			// instead: 
			// https://stackoverflow.com/questions/682438/hash-function-providing-unique-uint-from-an-integer-coordinate-pair
			// https://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
			// and arbitrarily ignore directionality
			std::size_t pairing = .5 * (x.first.id + x.second.id) * (x.first.id + x.second.id + 1) + x.second.id;
			return std::hash<std::size_t>()(pairing);
		}
	};

class Node{
	public:
		Node();
		Node(int id, std::string seq);
		std::string node_seq;
		int node_id;
		int bubble_id;
		int chain_id;
		int allele;
		bool visited;
		std::vector<std::pair<int, int>> childrenleft;
		std::vector<std::pair<int, int>> childrenright;
		bool operator==(const Node& other) const;
		bool operator!=(const Node& other) const;
		bool operator<(const Node& other) const;
		bool operator>(const Node& other) const;

};

class Bubble{
	public:
		Bubble();
		Bubble(Node source, Node sink, std::vector<Node> innerNodes);
		int id;
		Node source;
		Node sink;
		std::vector<Node> innerNodes;
		std::vector<Node> getNodes();
		std::vector<int> getNodeIds();
		Node getNode(int node_id);
};

class Chain{
	public:
		Chain();
		int id;
		std::vector<Bubble> bubbles;
		void addBubble(Bubble b);
		std::vector<Node> listChain();
		int size();
		
};

class Graph {
public:
	Graph();
	static Graph ReadGraph(std::string filename);
	void addNode(Node node);
	Node getNode(int node_id);
	void findBubbles();	
	void findBubble(Node node, bool direction, Chain* bchain);
	std::unordered_map<int, Node> nodes;
	std::unordered_map<DirectedNode, std::vector<DirectedNode>> edges;
	std::pair<DirectedNode, DirectedNode> getEdge(Node first, Node second);
	std::unordered_map<std::pair<DirectedNode, DirectedNode>, size_t> offsets;
	std::vector<Chain> chains;
	Chain getChain(int chain_id);
	void printChain(int chain_id);

};

#endif
