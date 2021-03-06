#ifndef ALIGNMENTREADER_HPP
#define ALIGNMENTREADER_HPP

#include <string>
#include <vector>
#include <algorithm>
#include <istream>
#include <ostream>
#include <unordered_map>
#include "graph.hpp"


/** 
* Reads a graph alignment file
**/

class AlignmentPath{
	public:
		AlignmentPath();
		AlignmentPath(std::vector<std::string> nodes);
		std::vector<std::string> nodes;
		int startpos;
		int endpos;
		float id;
		std::string name;
		int size();
		void print();
		std::vector<int> getRawIds();
		int firstnodelength;

};


class AlignmentReader{
	public:
		AlignmentReader();
		void readAlignmentfile(std::string filename, Graph graph);
		std::unordered_map<int, std::vector<AlignmentPath>> alignments;
};

#endif