#include "alignmentreader.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <vector>
#include <utility>
#include <set> 
#include <unordered_set> 
#include <regex>


using namespace std;


AlignmentPath::AlignmentPath():
nodes(),
startpos(0),
endpos(0),
id(0),
name(),
firstnodelength(0)
{
}

AlignmentPath::AlignmentPath(std::vector<std::string> nodes):
nodes(nodes),
startpos(0),
endpos(0),
id(0),
name(),
firstnodelength(0)
{
}

int AlignmentPath::size(){
	return nodes.size();
}

void AlignmentPath::print() {
	cout << "startpos: " << startpos << ", endpos: " << endpos << endl;
	for (auto& node: nodes){
		cout << node << ',';	
	}
	cout << endl;
}

//TODO put this function central somewhere
int raw_node_id(string s) {
//remove all non-numeric characters from s
s.erase(remove_if(s.begin(), s.end(), [](char c) { return !isdigit(c); }), s.end());
return stoi(s);
}

vector<int> AlignmentPath::getRawIds() {
	vector<int> result;
	for (auto& node: nodes) {
		result.push_back(raw_node_id(node));	
	}
	return(result);
}

AlignmentReader::AlignmentReader():
alignments()
{
}

void AlignmentReader::readAlignmentfile(std::string filename, Graph graph) {
	ifstream file(filename);
	string line;
	int counter = 0;
	ofstream myfile;
	string identityfile = filename.substr(0,filename.find(".gaf"));
	myfile.open (identityfile+"-alignment_identities.txt");

		
	while (getline(file, line)){
		
		if (!file.good()) break;
		if (line.size() == 0) continue;
		counter += 1;
	
		stringstream sstr {line};
		string name;
		string d1;
		string path;
		sstr >> name;
		//TODO dummy for HG00733_ontul_to_dip.r.utg.gaf
	//	string dum;
	//	sstr >> dum;
		
		sstr >> d1;
		sstr >> d1;
		sstr >> d1;
		sstr >> d1;
		sstr >> path;
		string dummy;
		sstr >> dummy;
		string startpos;
		string endpos;
		sstr >> startpos;
		sstr >> endpos;
		string d2;
		sstr >> d2;
		string length;
		sstr >> length;
	//	sstr >> d2;
		sstr >> d2;
		sstr >> d2;
		sstr >> d2;
		string id_string;
		//TODO dummy for mbg graph alignment
		string mbgdummy;
		sstr >> mbgdummy;
		sstr >> id_string;
		
		string id_delimiter = ":";
		size_t posi = 0;
		std::string token;
		vector<string> parts;
		while ((posi = id_string.find(id_delimiter)) != std::string::npos) {
		    token = id_string.substr(0, posi);
    		parts.push_back(token);
    id_string.erase(0, posi + id_delimiter.length());
}
		parts.push_back(id_string);
		assert(parts.at(0) == "id");

	//	string first_id = id_string.substr(0, id_string.find(id_delimiter));
	//	assert(first_id == "id");
	//	string id_value = id_string.substr(2, id_string.find(id_delimiter));
		string id_value = id_string;	
//		cout << "id_value: " << id_value << endl;
		float id_val = stof(id_value);

		
		vector<string> nodes;
		
		string const delims{ "<>" };

		size_t beg, pos = 0;
		string dir = "";
		vector<string> directions;
		while ((beg = path.find_first_not_of(delims, pos)) != string::npos){
			pos = path.find_first_of(delims, beg + 1);
			dir = path.substr(beg-1,1);
			directions.push_back(dir);
			nodes.push_back(path.substr(beg, pos - beg));
		}	
		myfile << name << "\t" << id_val << "\t";
		for (auto& node: nodes) {
			myfile << node << ',';		
		}
		myfile << "\t" << length;
		myfile << endl;
		assert(nodes.size() == directions.size());

	AlignmentPath alpath(nodes);		
		
		alpath.firstnodelength = -1;
		if (nodes.size() > 1) {
			DirectedNode from_firstnode = DirectedNode(raw_node_id(nodes[0]), directions[0] == ">");
			DirectedNode to_secondnode = DirectedNode(raw_node_id(nodes[1]), directions[1] == ">");
			int overlap = graph.offsets[make_pair(from_firstnode, to_secondnode)];
		//	TODO: add exception if node is not in graph		
		//	int firstnode_uniquelength = graph.getNode(raw_node_id(nodes[0])).node_seq.size() - overlap;
		int firstnode_uniquelength = 0;
			alpath.firstnodelength = firstnode_uniquelength;
		}
		
		alpath.startpos = stoi(startpos);
		alpath.endpos = stoi(endpos);
		alpath.id = id_val;
		alpath.name = name;
		for (auto& node: nodes) {
			if (node.length() > 0 ) {
				
			int nodeid = raw_node_id(node);
			int chain = graph.nodes[nodeid].chain_id;
			alignments[chain].push_back(alpath);
			}
		}
	}
	myfile.close();
}

/*
	path = line.strip().split('\t')[5]
	nodes = re.split('<|>', path)
	
#	clean_nodes = [int(node.replace("utg","").replace("l","")) for node in nodes if len(node)>= 1]
	clean_nodes = [int(re.sub('[^0-9]','', node)) for node in nodes if len(node)>= 1]
	node_id = clean_nodes[0]
	
	bubbleset = set()	
	chainset = set()
	alignedpath = AlignmentPath(clean_nodes)
	startposition = int(line.strip().split('\t')[7])
	endposition = int(line.strip().split('\t')[8])
	alignedpath.startpos = startposition
	alignedpath.endpos = endposition
	align.add_path(alignedpath)
	start_chain = graph.get_node(str(clean_nodes[0])).chain_id
	
	if (start_chain != False):
		start_align_dict[start_chain].append(alignedpath)
	for cleannode in clean_nodes:
		c = graph.get_node(str(cleannode)).chain_id
		b = graph.get_node(str(cleannode)).bubble_id
		if (c != False):
			align_dict[c].add(alignedpath)
		bubbleset.add(b)
		chainset.add(c)
	chainsets.append(chainset)
	bubblesets.append(bubbleset)	
*/

