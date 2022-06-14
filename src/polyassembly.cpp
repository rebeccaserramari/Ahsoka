#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "argumentparser.hpp"

#include "graph.hpp"
#include <jellyfish/mer_dna.hpp>
#include "/home/rebecca/work/whatshap-code/whatshap/src/read.h"
#include "alignmentreader.hpp"
#include "chainstoreadset.cpp"
#include "alignmentstoreadset.cpp"
//#include "jellyfishkmercounter.hpp"
#include <bits/stdc++.h>
//#include "stdlogger.h"
#include <thread>
#include <mutex>

using namespace std; 

void accumulator_function2(AlignmentReader& alignmentreader, Graph& graph, unordered_map<int, unordered_map<int, vector<vector<int>>>>& pathToAlleles, string readsetfile, bool shell_logging, vector<pair<int, int>>& size_sorting);
void sum_funcs(AlignmentReader& alignmentreader, Graph& graph, unordered_map<int, unordered_map<int, vector<vector<int>>>>& pathToAlleles, string readsetfile, bool shell_logging, vector<pair<int, int>>& size_sorting);                            

int main(int argc, char* argv[]) {

//	Logger logger = make_shared<StdLogger>();
	string gfafile = "";
	string alignmentfile = "";
	string readfile = "";
	string shortreadfile = "";
	string allelefile = "";
	string kmerfile = "";
	string countfile = "";
	string outfolder = "";
	int threads = 1;

	cerr << "Ahsoka: Haplotype assembly for diploid and polyploid genomes based on HiFi and ultra-long ONT data" << endl;
	cerr << "author: Rebecca Serra Mari" << endl;

	ArgumentParser argparser;
	argparser.add_command("Ahsoka [options] -g <graph.gfa> -a <alignments.gaf>");

	argparser.add_subcommand("only-bubbles", {'g', 'o'},{'t'});
	argparser.add_subcommand("phase", {'g', 'a', 'o'}, {'s','t'});
	
	string cmd;
	cmd = argparser.get_subcommand(argc, argv);
	if (cmd == "phase") {
		argparser.add_mandatory_argument('g', "genome assembly graph in gfa format, e.g. by hifiasm");
		argparser.add_mandatory_argument('a', "alignments of ONT reads to the assembly graph, in gaf format");
		argparser.add_mandatory_argument('o', "output folder to store output files");
		argparser.add_optional_argument('s',"", "additional long range phasing information (StrandSeq)");	
		argparser.add_optional_argument('t',"1", "number of threads to use");	
	}
	else if (cmd == "only-bubbles") {
		argparser.add_mandatory_argument('g', "genome assembly graph in gfa format, e.g. by hifiasm");		
		argparser.add_mandatory_argument('o', "output folder");	
		argparser.add_optional_argument('t',"1", "number of threads to use");	
	}

//argparser.add_optional_argument('r', "", "file with graph unitig sequences in fastq or fasta format");
//	argparser.add_optional_argument('s', "", "file with short read data, in fastq or fasta format");
//	//only for testing
	//argparser.add_optional_argument('t', "", "file with the allelepathsequence in fasta format");
	argparser.add_optional_argument('k', "", "kmerfile to be written to during the kmer counting");
	argparser.add_optional_argument('c', "", "outfile for summed up unique kmer counts in short read sample");
	try {
			argparser.parse(argc, argv);
	}
	catch (const runtime_error& e) {
		argparser.print_help();
		cerr << e.what() << endl;
		return 1;
	}
	catch (const exception& e) {
		return 0;
	}

	gfafile = argparser.get_arg_parameter('g');
	alignmentfile = argparser.get_arg_parameter('a');
	outfolder = argparser.get_arg_parameter('o');
	threads = stoi(argparser.get_arg_parameter('t'));
	//read in graph
  	Graph graph = Graph::ReadGraph(gfafile);
  	
	cout << "number of threads used: " << threads << endl;
	cout << "threads available: " << thread::hardware_concurrency() << endl;
	cout << "Step 1: Graph with " << graph.nodes.size() << " nodes read" << endl;
	
	//fill bubble information
	graph.findBubbles();	
	
	cout << "Step 2: Bubbles read" << endl;
	cout << "Number of bubble chains: " << graph.chains.size() << endl;	
	
	//outfile: take stem of graph file name
	ofstream bubblefile;
//	string infofile = gfafile.substr(0,gfafile.find(".gfa"));    
	string infofile = outfolder;
	bubblefile.open (infofile+"-bubbleinfo.txt");
	for (auto& chain: graph.chains) {
		bubblefile << "chain id: " << chain.id << "size: " << graph.getChain(chain.id).bubbles.size() << endl;
		for (auto& bubble: graph.getChain(chain.id).bubbles) {
			bubblefile << "bubble id: " << bubble.id << endl;
			bubblefile << "node id: " ;
			for (auto& node: bubble.getNodes())	{
				bubblefile << node.node_id << ",";
			}
			bubblefile << endl;
		}	
  	}
	bubblefile.close();
	//with option only-bubbles set, the program stops after writing the bubble file	
	if (cmd == "only-bubbles") {
		return 0;	
	}

	AlignmentReader alignmentreader;
	alignmentreader.readAlignmentfile(alignmentfile, graph);
	
	cout << "Step 3: Alignments read" << endl;	
	cout << "Number of alignments: " << alignmentreader.alignments.size() << endl;	
	//logger->log_info("Number of alignments: " + to_string(alignmentreader.alignments.size()));	

	//unordered_map<int, vector<vector<int>>> pathToAlleles;
	//pathToAlleles = ChainsToReadset(graph);	
	
	unordered_map<int, unordered_map<int,vector<vector<int>>>> chainpathToAlleles;
	chainpathToAlleles = ChainsToReadsetDetailed(graph);
	cout << "Step 4: Chain paths computed " << endl;
	cout << "Number of chain paths: " << chainpathToAlleles.size() <<  endl;
	
//	string readsetfile = alignmentfile.substr(0,alignmentfile.find(".gaf"));
	string readsetfile = outfolder;
	bool logging = false;
	
	//sort pathToAlleles by size of bubble chain, start with the longest ones first
	vector<pair<int, int>> size_sorting;
	for (auto& chainmap : chainpathToAlleles) {
		size_sorting.emplace_back(chainmap.second.size(), chainmap.first);
	}
	sort(begin(size_sorting), end(size_sorting), greater<>());
	
/*	vector<pair<int,int>> split_size_sorting_1;
	vector<pair<int,int>> split_size_sorting_2;
	copy(size_sorting.begin(), size_sorting.begin()+size_sorting.size()/2, back_inserter(split_size_sorting_1));
	copy(size_sorting.begin()+size_sorting.size()/2, size_sorting.end(), back_inserter(split_size_sorting_2));
	cout << "size sorting: " << size_sorting.size() << endl;
	cout << "split size sorting 1: " << split_size_sorting_1.size() << endl;
	cout << "split size sorting 2: " << split_size_sorting_2.size() << endl;
	
	for (int thr = 0; thr < threads; thr ++) {
		if (thr == 0) {
			thread t(alignmentsToReadset, ref(graph), ref(chainpathToAlleles), readsetfile, logging, ref(split_size_sorting_1));
			t.join();
		}
		else {
			thread t(alignmentsToReadset, ref(graph), ref(chainpathToAlleles), readsetfile, logging, ref(split_size_sorting_2));
			t.join();
		}	
	}*/
	//alignmentsToReadset(alignmentreader, graph, chainpathToAlleles, readsetfile,logging, size_sorting);	
	cout << "Step 5: Phasing processed" << endl;
	
	if (threads > 1) {
		cout << "multithread version" << endl;
	sum_funcs(alignmentreader, graph, chainpathToAlleles, readsetfile,logging, size_sorting);
	}
	else {
		cout << "single thread" << endl;
		cout << "size sorting: " << size_sorting.size() << endl;
		std::mutex g_display_mutex;
		alignmentsToReadset(alignmentreader, graph, chainpathToAlleles, readsetfile,logging, size_sorting,g_display_mutex);
		}	


	return 0;
}

void accumulator_function2(AlignmentReader& alignmentreader, Graph& graph, unordered_map<int, unordered_map<int, vector<vector<int>>>>& pathToAlleles, string readsetfile, bool shell_logging, vector<pair<int, int>>& size_sorting)
{
    /*acm = 0;
    for (unsigned int i = beginIndex; i < endIndex; ++i)
    {
        acm += v[i];
    }*/
    std::mutex g_display_mutex;
	//std::lock_guard<std::mutex> guard(g_display_mutex);
    alignmentsToReadset(alignmentreader, graph, pathToAlleles, readsetfile, shell_logging, size_sorting, g_display_mutex);	
}

void sum_funcs(AlignmentReader& alignmentreader, Graph& graph, unordered_map<int, unordered_map<int, vector<vector<int>>>>& pathToAlleles, string readsetfile, bool shell_logging, vector<pair<int, int>>& size_sorting){
			vector<pair<int,int>> split_size_sorting_1;
		//	copy(size_sorting.begin(), size_sorting.begin()+size_sorting.size()/2, back_inserter(split_size_sorting_1));
			split_size_sorting_1.push_back(size_sorting.at(0));
			split_size_sorting_1.push_back(size_sorting.at(2));
			split_size_sorting_1.push_back(size_sorting.at(4));
			split_size_sorting_1.push_back(size_sorting.at(6));
			split_size_sorting_1.push_back(size_sorting.at(8));
			
			vector<pair<int,int>> split_size_sorting_2;
		//	copy(size_sorting.begin()+size_sorting.size()/2, size_sorting.end(), back_inserter(split_size_sorting_2));
			split_size_sorting_2.push_back(size_sorting.at(1));
			split_size_sorting_2.push_back(size_sorting.at(3));
			split_size_sorting_2.push_back(size_sorting.at(5));
			split_size_sorting_2.push_back(size_sorting.at(7));
			split_size_sorting_2.push_back(size_sorting.at(9));
			
			cout << "size sorting: " << size_sorting.size() << endl;
			cout << "split size sorting 1: " << split_size_sorting_1.size() << endl;
			cout << "split size sorting 2: " << split_size_sorting_2.size() << endl;
			   
        	std::thread t1(accumulator_function2, ref(alignmentreader), ref(graph), ref(pathToAlleles),readsetfile, shell_logging,ref(split_size_sorting_1));
        	std::thread t2(accumulator_function2, ref(alignmentreader), ref(graph), ref(pathToAlleles),readsetfile, shell_logging,ref(split_size_sorting_2));
     			
        t1.join();
        t2.join();
			
    //    std::cout << "acm1: " << acm1 << endl;
    //    std::cout << "acm2: " << acm2 << endl;
    //    std::cout << "acm1 + acm2: " << acm1 + acm2 << endl;

 
 }