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

using namespace std; 

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

	cerr << "Ahsoka: Haplotype assembly for diploid and polyploid genomes based on HiFi and ultra-long ONT data" << endl;
	cerr << "author: Rebecca Serra Mari" << endl;

	ArgumentParser argparser;
	argparser.add_command("Ahsoka [options] -g <graph.gfa> -a <alignments.gaf>");

	argparser.add_subcommand("only-bubbles", {'g', 'o'},{});
	argparser.add_subcommand("phase", {'g', 'a', 'o'}, {'s'});
	
	string cmd;
	cmd = argparser.get_subcommand(argc, argv);
	if (cmd == "phase") {
		argparser.add_mandatory_argument('g', "genome assembly graph in gfa format, e.g. by hifiasm");
		argparser.add_mandatory_argument('a', "alignments of ONT reads to the assembly graph, in gaf format");
		argparser.add_mandatory_argument('o', "output folder to store output files");
		argparser.add_optional_argument('s',"", "additional long range phasing information (StrandSeq)");	
	}
	else if (cmd == "only-bubbles") {
		argparser.add_mandatory_argument('g', "genome assembly graph in gfa format, e.g. by hifiasm");		
		argparser.add_mandatory_argument('o', "output folder");		
	}

//argparser.add_optional_argument('r', "", "file with graph unitig sequences in fastq or fasta format");
//	argparser.add_optional_argument('s', "", "file with short read data, in fastq or fasta format");
//	//only for testing
	argparser.add_optional_argument('t', "", "file with the allelepathsequence in fasta format");
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
	
	//read in graph
  	Graph graph = Graph::ReadGraph(gfafile);
	
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
	alignmentsToReadset(alignmentreader, graph, chainpathToAlleles, readsetfile,logging);	
	cout << "Step 5: Phasing processed" << endl;
	


	return 0;
}



