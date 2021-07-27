#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "argumentparser.hpp"
#include "timer.hpp"

#include "graph.hpp"
#include <jellyfish/mer_dna.hpp>
#include "/home/rebecca/work/whatshap-code/whatshap/src/read.h"
#include "alignmentreader.hpp"
#include "chainstoreadset.cpp"
#include "alignmentstoreadset.cpp"
#include "jellyfishkmercounter.hpp"
#include "uniquekmers.hpp"
#include <bits/stdc++.h>

using namespace std;


int EditDistDP(string str1, string str2)
{
    int len1 = str1.length();
    int len2 = str2.length();
 
    // Create a DP array to memoize result
    // of previous computations
    int DP[2][len1 + 1];
 
    // To fill the DP array with 0
    memset(DP, 0, sizeof DP);
 
    // Base condition when second string
    // is empty then we remove all characters
    for (int i = 0; i <= len1; i++)
        DP[0][i] = i;
 
    // Start filling the DP
    // This loop run for every
    // character in second string
    for (int i = 1; i <= len2; i++) {
        // This loop compares the char from
        // second string with first string
        // characters
        for (int j = 0; j <= len1; j++) {
            // if first string is empty then
            // we have to perform add character
            // operation to get second string
            if (j == 0)
                DP[i % 2][j] = i;
 
            // if character from both string
            // is same then we do not perform any
            // operation . here i % 2 is for bound
            // the row number.
            else if (str1[j - 1] == str2[i - 1]) {
                DP[i % 2][j] = DP[(i - 1) % 2][j - 1];
            }
 
            // if character from both string is
            // not same then we take the minimum
            // from three specified operation
            else {
                DP[i % 2][j] = 1 + min(DP[(i - 1) % 2][j],
                                       min(DP[i % 2][j - 1],
                                           DP[(i - 1) % 2][j - 1]));
            }
        }
    }
 
    // after complete fill the DP array
    // if the len2 is even then we end
    // up in the 0th row else we end up
    // in the 1th row so we take len2 % 2
    // to get row
    cout << DP[len2 % 2][len1] << endl;
    
    return(DP[len2 % 2][len1]);
}
 

int my_raw_id(string s) {
//remove all non-numeric characters from s
s.erase(remove_if(s.begin(), s.end(), [](char c) { return !isdigit(c); }), s.end());
return stoi(s);
}

int main(int argc, char* argv[]) {

	Timer timer;
	double time_read_graph;
	string gfafile = "";
	string alignmentfile = "";
	string readfile = "";
	string shortreadfile = "";
	string allelefile = "";
	string kmerfile = "";
	string countfile = "";

	cerr << "Phoenix: Polyploid genome assembly based on HiFi, ONT and short-read data" << endl;
	cerr << "Rebecca Serra Mari" << endl;

	ArgumentParser argparser;
	argparser.add_command("Polyassembly [options] -g <graph.gfa> -a <alignments.gaf>");

	argparser.add_subcommand("only-bubbles", {'g'},{});
	argparser.add_subcommand("phase", {'g', 'a'}, {'s'});
	
	//argparser.allow_argument('g', "eval");
	//argparser.allow_argument('g', "only-bubbles");
	//argparser.allow_argument('g', "phase");

	string cmd;
	cmd = argparser.get_subcommand(argc, argv);
	if (cmd == "phase") {
		argparser.add_mandatory_argument('g', "genome assembly graph in gfa format, e.g. by hifiasm");
		argparser.add_mandatory_argument('a', "alignments of ONT reads to the assembly graph, in gaf format");
		argparser.add_optional_argument('s',"", "additional long range phasing information (StrandSeq)");	
	}
	else if (cmd == "only-bubbles") {
		argparser.add_mandatory_argument('g', "genome assembly graph in gfa format, e.g. by hifiasm");		
	}

//argparser.add_optional_argument('r', "", "file with graph unitig sequences in fastq or fasta format");
//	argparser.add_optional_argument('s', "", "file with short read data, in fastq or fasta format");
//	//only for testing
	argparser.add_optional_argument('t', "", "file with the allelepathsequence in fasta format");
	argparser.add_optional_argument('k', "", "kmerfile to be written to during the kmer counting");
	argparser.add_optional_argument('c', "", "outfile for summed up unique kmer counts in short read sample");
	cout << "before argparse " << endl;
		try {
			argparser.parse(argc, argv);
		} catch (const runtime_error& e) {
			argparser.print_help();
			cerr << e.what() << endl;
			return 1;
		} catch (const exception& e) {
			return 0;
		}

	gfafile = argparser.get_arg_parameter('g');
	alignmentfile = argparser.get_arg_parameter('a');
	cout << "gfafile: " << gfafile << endl;
//	readfile = argparser.get_argument('r');
//	shortreadfile = argparser.get_argument('s');
//	allelefile = argparser.get_argument('t');
//	kmerfile = argparser.get_argument('k');
//	countfile = argparser.get_argument('c');

	//cerr << "Files and parameters used:" << endl;
	//argparser.info();

/*	argparser.add_subcommand("eval");
	argparser.add_subcommand("only-bubbles");
	argparser.add_subcommand("phase");
	argparser.allow_argument('g', "eval");
	argparser.allow_argument('g', "only-bubbles");
	argparser.allow_argument('g', "phase");
	

	string cmd;
	cmd = argparser.get_subcommand(argv);
	//TODO check whether the subcommand list can be added when filling the add_argument functions: in that function, add the arg to the allowed commands
	//then, in the check below, only the function call is necessary?
	//does not account for args that are mandatory in one case, but optional in another!
	//add for each command individually: argparser.add_mandatory_argument("eval", 'g', "genome assembly graph in gfa format, e.g. by hifiasm"); argparser.add_optional_argument("something", 'g', "genome assembly graph in gfa format, e.g. by hifiasm");
	if (cmd == "eval") {
		//start evaluation function
		//eval_phasing(phasedfile);
		string resultfile = "";
		string gfafile = "";
		argparser.add_mandatory_argument('r',"result file with phased node sequences and switch distances");
		argparser.add_mandatory_argument('g',"genome assembly graph in gfa format, e.g. by hifiasm");
				try {
		argparser.parse(argc, argv);
	} catch (const runtime_error& e) {
		argparser.usage();
		cerr << e.what() << endl;
		return 1;
	} catch (const exception& e) {
		return 0;
	}

		resultfile = argparser.get_argument('r');
		gfafile = argparser.get_argument('f');
		cout << "gfafile: " << gfafile << endl;
	}
	else if (cmd == "only-bubbles") {
		string gfafile = "";
		string outfile = "";
		argparser.add_mandatory_argument('g',"genome assembly graph in gfa format, e.g. by hifiasm");
		argparser.add_optional_argument('o',"cout","output file for the bubbleinfo, otherwise cout");
				try {
		argparser.parse(argc, argv);
	} catch (const runtime_error& e) {
		argparser.usage();
		cerr << e.what() << endl;
		return 1;
	} catch (const exception& e) {
		return 0;
	}

		gfafile = argparser.get_argument('g');
		outfile = argparser.get_argument('o');
		cout << "gfafile b: " << gfafile << endl;
		//function for writing the bubble file	
	}
	else if (cmd == "phase") {
		string gfafile = "";
		string outfile = "";
		argparser.add_mandatory_argument('g',"genome assembly graph in gfa format, e.g. by hifiasm");
		argparser.add_mandatory_argument('a',"alignments of ONT reads to the assembly graph, in gaf format");
		argparser.add_optional_argument('o',"cout","output file for the bubbleinfo, otherwise cout");
				try {
		argparser.parse(argc, argv);
	} catch (const runtime_error& e) {
		argparser.usage();
		cerr << e.what() << endl;
		return 1;
	} catch (const exception& e) {
		return 0;
	}

		gfafile = argparser.get_argument('g');
		alignmentfile = argparser.get_argument('a');
		outfile = argparser.get_argument('o');
		cout << "gfafile b: " << gfafile << endl;
		//function for writing the bubble file	
	}
*/	

  	Graph graph = Graph::ReadGraph(gfafile);
	
	cout << "Graph read " << endl;
	cout << graph.nodes.size() << endl;

	
	graph.findBubbles();	
	
	cout << "bubbles read " << endl;
	cout <<"no. of nodes: " << graph.nodes.size() << endl;
	cout << "number of chains: " << graph.chains.size() << endl;	

/*
//hifiasm normal	
//	vector<string> hap1 = {"utg000288l", "utg001036l", "utg000554l", "utg001197l", "utg000519l", "utg001171l", "utg000309l", "utg001111l", "utg000577l", "utg001296l", "utg000561l", "utg000148l", "utg001295l", "utg000118l", "utg000516l", "utg001138l", "utg000127l", "utg000446l", "utg000942l", "utg000954l", "utg000564l", "utg001299l", "utg000124l", "utg001144l", "utg000575l", "utg001157l", "utg000176l", "utg001177l", "utg000763l", "utg000487l", "utg001112l", "utg000533l", "utg000387l", "utg001347l", "utg000424l", "utg001381l", "utg000815l", "utg001285l", "utg000336l", "utg001091l", "utg000167l", "utg001236l", "utg000534l", "utg001018l", "utg000339l", "utg001185l", "utg000088l", "utg001045l", "utg000679l", "utg001094l", "utg000262l", "utg000332l", "utg000342l", "utg000999l", "utg000278l", "utg001311l", "utg000308l", "utg001211l", "utg000155l", "utg000492l", "utg000571l", "utg001272l", "utg000888l", "utg001101l", "utg000185l", "utg000078l", "utg001099l", "utg000284l", "utg000084l", "utg001267l", "utg000670l", "utg000749l", "utg000018l", "utg000298l", "utg001255l", "utg000592l", "utg000905l", "utg000195l", "utg001006l", "utg000279l", "utg001067l", "utg000196l", "utg001415l", "utg000483l", "utg000179l", "utg000276l", "utg001389l", "utg000199l", "utg000355l", "utg000171l", "utg000347l", "utg001051l", "utg001151l", "utg000648l", "utg000477l", "utg001121l", "utg000719l", "utg000744l", "utg000574l", "utg001254l", "utg000596l", "utg001360l", "utg000740l", "utg001393l", "utg001382l", "utg001206l", "utg001086l", "utg000716l", "utg000632l", "utg001262l", "utg000526l", "utg000780l", "utg001060l", "utg001335l", "utg001148l", "utg000524l", "utg001351l", "utg000103l", "utg001434l", "utg000542l", "utg001145l", "utg000568l", "utg000878l", "utg000068l", "utg001373l", "utg000767l", "utg000286l", "utg000165l", "utg000700l", "utg001076l", "utg000220l", "utg000992l", "utg001048l", "utg001245l", "utg000981l", "utg001069l", "utg000718l", "utg000727l", "utg000984l", "utg000313l", "utg001141l", "utg000741l", "utg001325l", "utg000197l", "utg001013l", "utg000739l", "utg000747l", "utg000929l", "utg000566l", "utg001203l", "utg000772l", "utg001092l", "utg000786l", "utg000931l", "utg000952l", "utg000985l", "utg000488l", "utg001059l", "utg000495l", "utg001340l", "utg000544l", "utg000664l", "utg001167l", "utg000484l", "utg001188l", "utg000192l", "utg001029l", "utg000113l", "utg000345l", "utg000439l", "utg000093l", "utg000683l", "utg000944l", "utg001127l", "utg001320l", "utg000045l", "utg001187l", "utg000301l", "utg000426l", "utg000249l", "utg001195l", "utg000462l", "utg000095l", "utg001223l", "utg000131l", "utg000143l", "utg001330l", "utg000010l", "utg000188l", "utg000375l", "utg000831l", "utg001279l", "utg000177l", "utg000187l", "utg001058l", "utg000651l", "utg001140l", "utg000049l", "utg001406l", "utg000482l", "utg000476l", "utg001288l", "utg000472l", "utg001234l", "utg000491l", "utg001186l", "utg000264l", "utg001242l", "utg000365l", "utg001413l", "utg000117l", "utg000303l", "utg000052l", "utg000833l", "utg000855l", "utg000766l", "utg001082l", "utg000330l", "utg001341l", "utg000692l", "utg001361l", "utg001095l", "utg001213l", "utg000151l", "utg001166l", "utg000920l", "utg000918l", "utg000134l", "utg001392l", "utg000557l", "utg000076l", "utg001378l", "utg000530l", "utg001260l", "utg000273l", "utg001194l", "utg000075l", "utg000922l", "utg001230l", "utg000960l", "utg001326l", "utg000994l", "utg000986l", "utg000503l", "utg001225l", "utg001072l", "utg000297l", "utg001070l", "utg000759l", "utg001283l", "utg000928l", "utg001377l", "utg000133l", "utg000428l", "utg000479l", "utg001253l", "utg000361l", "utg001269l", "utg000770l", "utg001327l", "utg000475l", "utg001315l", "utg000714l", "utg001204l", "utg000694l", "utg001041l", "utg000377l", "utg000444l", "utg000417l", "utg001074l", "utg000724l", "utg001039l", "utg000569l", "utg000011l", "utg001286l", "utg000474l", "utg001324l", "utg000385l", "utg001243l", "utg000217l", "utg000639l", "utg001431l", "utg000972l", "utg000968l", "utg000976l", "utg001020l", "utg000966l", "utg000995l", "utg000209l", "utg001418l", "utg000836l", "utg001224l", "utg000392l", "utg001201l", "utg001049l", "utg001117l", "utg000742l", "utg001043l", "utg001085l", "utg000567l", "utg000379l", "utg001416l", "utg000455l", "utg001400l", "utg000473l", "utg000496l", "utg000587l", "utg001210l", "utg000147l", "utg001159l", "utg000182l", "utg001372l", "utg000162l", "utg001089l", "utg000642l", "utg001142l", "utg000646l", "utg000549l", "utg000854l", "utg001412l", "utg000316l", "utg001083l", "utg000609l", "utg001080l", "utg000354l", "utg001077l", "utg000586l", "utg000430l", "utg000512l", "utg000126l", "utg000404l", "utg000490l", "utg000835l", "utg000853l", "utg001343l", "utg000774l", "utg001155l", "utg000322l", "utg000431l", "utg000140l", "utg000407l", "utg000122l", "utg000969l"};
//	vector<string> hap2 = {"utg000493l", "utg001036l", "utg000108l", "utg001197l", "utg000004l", "utg001171l", "utg000800l", "utg001111l", "utg000200l", "utg001296l", "utg000058l", "utg000757l", "utg001295l", "utg000516l", "utg001138l", "utg000127l", "utg001264l", "utg000954l", "utg000400l", "utg001299l", "utg000738l", "utg001144l", "utg000302l", "utg001157l", "utg000612l", "utg001177l", "utg000321l", "utg000050l", "utg001112l", "utg000139l", "utg000627l", "utg001347l", "utg000818l", "utg000834l", "utg000815l", "utg001285l", "utg000604l", "utg001091l", "utg000559l", "utg001236l", "utg000071l", "utg001018l", "utg000704l", "utg001185l", "utg000540l", "utg001045l", "utg000232l", "utg001094l", "utg000610l", "utg000814l", "utg000816l", "utg000999l", "utg000680l", "utg001311l", "utg000550l", "utg001211l", "utg000753l", "utg000310l", "utg000274l", "utg000032l", "utg000023l", "utg001272l", "utg000983l", "utg001101l", "utg000812l", "utg000762l", "utg001099l", "utg000777l", "utg000525l", "utg001267l", "utg000237l", "utg000158l", "utg000486l", "utg000845l", "utg001255l", "utg000225l", "utg000905l", "utg000663l", "utg001006l", "utg000748l", "utg001067l", "utg000710l", "utg001415l", "utg000022l", "utg000655l", "utg000819l", "utg001389l", "utg000535l", "utg000631l", "utg001365l", "utg001151l", "utg000186l", "utg000107l", "utg001121l", "utg000211l", "utg000327l", "utg000381l", "utg000413l", "utg001254l", "utg000331l", "utg001360l", "utg000740l", "utg001393l", "utg001382l", "utg001206l", "utg001116l", "utg000317l", "utg000406l", "utg001262l", "utg000253l", "utg000201l", "utg001060l", "utg001252l", "utg001148l", "utg000072l", "utg001351l", "utg000811l", "utg001434l", "utg000112l", "utg001145l", "utg000003l", "utg000878l", "utg000690l", "utg001373l", "utg000388l", "utg000851l", "utg000460l", "utg000367l", "utg001076l", "utg000688l", "utg000992l", "utg001065l", "utg001245l", "utg001221l", "utg001069l", "utg000214l", "utg000224l", "utg000984l", "utg000599l", "utg001141l", "utg000312l", "utg001325l", "utg000597l", "utg001013l", "utg000401l", "utg000281l", "utg000929l", "utg000292l", "utg001203l", "utg000306l", "utg001092l", "utg000343l", "utg000931l", "utg000952l", "utg000985l", "utg000066l", "utg001059l", "utg000057l", "utg001340l", "utg000296l", "utg000397l", "utg001167l", "utg000006l", "utg001188l", "utg000708l", "utg001029l", "utg000603l", "utg000653l", "utg000093l", "utg000230l", "utg000944l", "utg001130l", "utg001320l", "utg000468l", "utg001187l", "utg000794l", "utg000809l", "utg000471l", "utg001195l", "utg000030l", "utg000465l", "utg001223l", "utg000517l", "utg000637l", "utg001330l", "utg000614l", "utg000605l", "utg000375l", "utg000831l", "utg001279l", "utg000728l", "utg000598l", "utg001058l", "utg000651l", "utg001140l", "utg000508l", "utg001406l", "utg000156l", "utg000348l", "utg001288l", "utg000102l", "utg001234l", "utg000146l", "utg001186l", "utg000454l", "utg001242l", "utg000640l", "utg001413l", "utg000478l", "utg000521l", "utg000501l", "utg000159l", "utg000271l", "utg001082l", "utg000685l", "utg001341l", "utg000041l", "utg001361l", "utg001095l", "utg001213l", "utg000705l", "utg001166l", "utg001216l", "utg000918l", "utg000615l", "utg001392l", "utg000132l", "utg000666l", "utg001378l", "utg000337l", "utg001260l", "utg000795l", "utg001194l", "utg000636l", "utg000922l", "utg001176l", "utg000960l", "utg001326l", "utg001106l", "utg000986l", "utg000269l", "utg001169l", "utg001072l", "utg000808l", "utg001070l", "utg000422l", "utg001283l", "utg001310l", "utg001377l", "utg000668l", "utg000428l", "utg000282l", "utg001253l", "utg000850l", "utg001269l", "utg000277l", "utg001327l", "utg000111l", "utg001315l", "utg000359l", "utg001204l", "utg000358l", "utg001041l", "utg000746l", "utg001409l", "utg001074l", "utg000724l", "utg001039l", "utg000378l", "utg000580l", "utg000678l", "utg000448l", "utg001286l", "utg000069l", "utg001324l", "utg000551l", "utg001243l", "utg000790l", "utg000036l", "utg001431l", "utg001035l", "utg001162l", "utg000976l", "utg000951l", "utg001322l", "utg000995l", "utg000513l", "utg001418l", "utg000323l", "utg001224l", "utg000677l", "utg001201l", "utg001278l", "utg001117l", "utg000017l", "utg001052l", "utg001085l", "utg000161l", "utg000720l", "utg001416l", "utg000166l", "utg001400l", "utg000025l", "utg000070l", "utg000136l", "utg001210l", "utg000578l", "utg001159l", "utg000527l", "utg001372l", "utg000660l", "utg001089l", "utg000215l", "utg001142l", "utg000037l", "utg000257l", "utg001412l", "utg000788l", "utg001083l", "utg000235l", "utg001080l", "utg000354l", "utg001077l", "utg000019l", "utg000430l", "utg000002l", "utg000412l", "utg000404l", "utg000395l", "utg000389l", "utg000490l", "utg000415l", "utg001343l", "utg000774l", "utg001155l", "utg000830l", "utg000847l", "utg000687l", "utg000778l", "utg000733l", "utg000969l"};
//	vector<int> hap1_switches = {1, 5, 7, 10, 13, 15, 16, 21, 23, 25, 26, 37, 41, 45, 47, 49, 58, 63, 70, 71, 74, 76, 82, 83, 92, 105, 106, 116, 118, 122, 135, 138, 144, 153, 155, 164, 170, 174, 180, 181, 188, 191, 194, 196, 198, 205, 212, 224, 228, 229, 231, 233, 243, 250, 255, 265, 269, 272, 274, 276, 279, 284, 295, 299, 301, 307, 310, 311, 313, 321, 325, 329, 330, 338, 340};
//	vector<int> hap2_switches = {1, 5, 7, 10, 14, 17, 21, 23, 25, 31, 39, 41, 45, 47, 56, 63, 70, 71, 74, 76, 82, 83, 88, 99, 105, 115, 117, 121, 126, 137, 139, 141, 143, 163, 164, 165, 167, 172, 178, 179, 183, 184, 196, 203, 210, 217, 225, 226, 228, 230, 240, 248, 254, 262, 267, 268, 272, 274, 277, 282, 287, 291, 293, 297, 299, 305, 311, 318, 319, 329, 330, 335}; 

//hifiasm trio
	vector<string> hap1 = {"utg000393l","utg002055l","utg000119l","utg002209l","utg000004l","utg002183l","utg000436l","utg002126l","utg000234l","utg002307l","utg000062l","utg002306l","utg001026l","utg002152l","utg000138l","utg001974l","utg000663l","utg002310l","utg000134l","utg002158l","utg000428l","utg002170l","utg000195l","utg002189l","utg000054l","utg002127l","utg000151l","utg000620l","utg002356l","utg000748l","utg002388l","utg001543l","utg002296l","utg000490l","utg002106l","utg000183l","utg002248l","utg000078l","utg002038l","utg000494l","utg002197l","utg000098l","utg002063l","utg000285l","utg002109l","utg000335l","utg000483l","utg000507l","utg002019l","utg000370l","utg002321l","utg000435l","utg002223l","utg000170l","utg000214l","utg000032l","utg000023l","utg002283l","utg001910l","utg002116l","utg000206l","utg000087l","utg002114l","utg000383l","utg000093l","utg002278l","utg001387l","utg001224l","utg000275l","utg000250l","utg000816l","utg000438l","utg000174l","utg000540l","utg000077l","utg000431l","utg000454l","utg000637l","utg000623l","utg000632l","utg000411l","utg000644l","utg000790l","utg000687l","utg000786l","utg000018l","utg000533l","utg000161l","utg000257l","utg000743l","utg000419l","utg000258l","utg000269l","utg000441l","utg000866l","utg000553l","utg000882l","utg000626l","utg000243l","utg000545l","utg000613l","utg000581l","utg000676l","utg000160l","utg000408l","utg000739l","utg000310l","utg000756l","utg000510l","utg000425l","utg000781l","utg000232l","utg000478l","utg000497l","utg000798l","utg000749l","utg000656l","utg000423l","utg000754l","utg002267l","utg000279l","utg001927l","utg000227l","utg002026l","utg000372l","utg002083l","utg000228l","utg002421l","utg000022l","utg002396l","utg000233l","utg000531l","utg000187l","utg002164l","utg000207l","utg000118l","utg002135l","utg000254l","utg000474l","utg000605l","utg000700l","utg002266l","utg000482l","utg002369l","utg001371l","utg002400l","utg002133l","utg002218l","utg000451l","utg000672l","utg002274l","utg000237l","utg002076l","utg002264l","utg002161l","utg000079l","utg002360l","utg000113l","utg002438l","utg000123l","utg002159l","utg000003l","utg001900l","utg000074l","utg002381l","utg000624l","utg000180l","utg000564l","utg002092l","utg000270l","utg002012l","utg002081l","utg002257l","utg002001l","utg002085l","utg000262l","utg000278l","utg002004l","utg000442l","utg002155l","utg000440l","utg002334l","utg000230l","utg002033l","utg000665l","utg000379l","utg001950l","utg000402l","utg002215l","utg001434l","utg002107l","utg000508l","utg001952l","utg001972l","utg002005l","utg000072l","utg002075l","utg000061l","utg002349l","utg000420l","utg000649l","utg002180l","utg000006l","utg002200l","utg000223l","utg002049l","utg000100l","utg001964l","utg002141l","utg002329l","utg000048l","utg002199l","utg000317l","utg002207l","utg000030l","utg002235l","utg000143l","utg000154l","utg002339l","utg000010l","utg000210l","utg000592l","utg001624l","utg002290l","utg000196l","utg000208l","utg002074l","utg001194l","utg002154l","utg000052l","utg002413l","utg000172l","utg000518l","utg002299l","utg000111l","utg002246l","utg000158l","utg002198l","utg000338l","utg002254l","utg000561l","utg002419l","utg000127l","utg000056l","utg000355l","utg002098l","utg000480l","utg002350l","utg000043l","utg002370l","utg002110l","utg002225l","utg000165l","utg002179l","utg001942l","utg001940l","utg000146l","utg002399l","utg000144l","utg000084l","utg002385l","utg000491l","utg002272l","utg000359l","utg002206l","utg000083l","utg001944l","utg002242l","utg001980l","utg002335l","utg002014l","utg002006l","utg000351l","utg002182l","utg002088l","utg000421l","utg002086l","utg000742l","utg002294l","utg000145l","utg000768l","utg000985l","utg002265l","utg000554l","utg002280l","utg000368l","utg002336l","utg000122l","utg002325l","utg000549l","utg002216l","utg000541l","utg002060l","utg000593l","utg000926l","utg000715l","utg002090l","utg001332l","utg002058l","utg000598l","utg000011l","utg002297l","utg000075l","utg002333l","utg000614l","utg002255l","utg000265l","utg000036l","utg002435l","utg002054l","utg002175l","utg001996l","utg002040l","utg001986l","utg002015l","utg000251l","utg002424l","utg000463l","utg002236l","utg000635l","utg002213l","utg002067l","utg002131l","utg000017l","utg002101l","utg000176l","utg002422l","utg000182l","utg002407l","utg000025l","utg000076l","utg000148l","utg002222l","utg000159l","utg002172l","utg000203l","utg002380l","utg000177l","utg002104l","utg000263l","utg002156l","utg000037l","utg000448l","utg002099l","utg000295l","utg002096l","utg000528l","utg002093l","utg000019l","utg000795l","utg000468l","utg000868l","utg000040l","utg001218l","utg000396l","utg000780l","utg000536l","utg000222l","utg000332l","utg000141l","utg000546l","utg000046l","utg000002l","utg000067l","utg000727l","utg000950l","utg000322l","utg000862l","utg000894l","utg000933l","utg000277l","utg000566l","utg000437l","utg000728l","utg000462l","utg000599l","utg000789l","utg000792l","utg000588l","utg000266l","utg000837l","utg000697l","utg000938l","utg000784l","utg000643l","utg000496l","utg000631l","utg000612l","utg000426l","utg000841l","utg000386l","utg000820l","utg000304l","utg000406l","utg000884l","utg000832l","utg000219l","utg000292l","utg000886l","utg000836l","utg000217l","utg000337l","utg000825l","utg000495l","utg000904l","utg000189l","utg000806l","utg000630l","utg000169l","utg000470l","utg000512l","utg000471l","utg000171l","utg000256l","utg000289l","utg000499l","utg000641l","utg000719l","utg000709l","utg000835l","utg002352l","utg001438l","utg002168l"};
	vector<string> hap2 = {};
	vector<int> hap1_switches = {11, 13, 15, 17, 30, 32, 142, 147, 220, 223, 226, 228, 249, 253, 282, 284, 296, 298, 335, 336, 343, 344, 352, 353};
	vector<int> hap2_switches = {};	
	
//me on the pgas graph
 //	vector<string> hap1 = {"utg000840l", "utg000088l", "utg001040l", "utg000135l", "utg000761l", "utg000018l", "utg000786l", "utg000730l", "utg000859l", "utg000288l", "utg000880l", "utg000247l", "utg001042l", "utg000174l", "utg000688l", "utg000664l", "utg000672l", "utg000718l", "utg001053l", "utg000036l", "utg000896l", "utg000287l", "utg000962l", "utg000066l", "utg000931l", "utg000721l", "utg000011l", "utg000253l", "utg000750l", "utg000130l", "utg000723l", "utg000514l", "utg000862l", "utg000271l", "utg000956l", "utg000334l", "utg000965l", "utg000562l", "utg000919l", "utg000273l", "utg000905l", "utg000339l", "utg000756l", "utg000494l", "utg001010l", "utg000636l", "utg000929l", "utg000556l", "utg000746l", "utg000231l", "utg000748l", "utg000681l", "utg000363l", "utg000831l", "utg000660l", "utg000964l", "utg000778l", "utg000885l", "utg000631l", "utg000071l", "utg000853l", "utg000216l", "utg000911l", "utg000256l", "utg001011l", "utg000411l", "utg001020l", "utg000461l", "utg000629l", "utg000873l", "utg000828l", "utg000520l", "utg000870l", "utg000769l", "utg000997l", "utg000512l", "utg000979l", "utg000504l", "utg000758l", "utg000338l", "utg001038l", "utg000479l", "utg000895l", "utg000313l", "utg000846l", "utg000352l", "utg000889l", "utg000331l", "utg000933l", "utg000336l", "utg001034l", "utg000368l", "utg000804l", "utg000485l", "utg000734l", "utg000449l", "utg000926l", "utg000594l", "utg000285l", "utg000454l", "utg000969l", "utg000376l", "utg000879l", "utg000854l", "utg000321l", "utg000324l", "utg000200l", "utg000847l", "utg000044l", "utg000960l", "utg000795l", "utg000650l", "utg000085l", "utg000713l", "utg000158l", "utg000848l", "utg000006l", "utg000829l", "utg000230l", "utg000978l", "utg000056l", "utg000735l", "utg000063l", "utg000680l", "utg000819l", "utg000639l", "utg000572l", "utg000766l", "utg000236l", "utg000861l", "utg000227l", "utg000637l", "utg000222l", "utg000700l", "utg000163l", "utg000963l", "utg000240l", "utg000805l", "utg000241l", "utg000679l", "utg000177l", "utg000745l", "utg000676l", "utg000898l", "utg000741l", "utg000686l", "utg000183l", "utg000752l", "utg000139l", "utg001006l", "utg000065l", "utg000607l", "utg000003l", "utg000809l", "utg000103l", "utg001055l", "utg000095l", "utg000988l", "utg000068l", "utg000812l", "utg000904l", "utg000736l", "utg000167l", "utg000913l", "utg000864l", "utg000244l", "utg000294l", "utg000788l", "utg001021l", "utg000107l", "utg000996l", "utg000252l", "utg000906l", "utg000175l", "utg000790l", "utg000098l", "utg000815l", "utg001017l", "utg000390l", "utg000732l", "utg000342l", "utg001039l", "utg000524l", "utg000743l", "utg000550l", "utg000697l", "utg000491l", "utg000624l", "utg000443l", "utg000907l", "utg000019l", "utg000917l", "utg000381l", "utg000773l", "utg000558l", "utg000775l", "utg000678l", "utg000922l", "utg000869l", "utg000354l", "utg000238l", "utg000952l", "utg000219l", "utg000691l", "utg000208l", "utg000768l", "utg000189l", "utg000726l", "utg000083l", "utg000845l", "utg000257l", "utg000704l", "utg000067l", "utg000891l", "utg000141l", "utg000765l", "utg000985l", "utg000300l", "utg000255l", "utg000120l", "utg000783l", "utg000049l", "utg000839l", "utg000148l", "utg000820l", "utg000234l", "utg000808l", "utg000109l", "utg000943l", "utg000293l", "utg000657l", "utg000802l", "utg000111l", "utg000915l", "utg000106l", "utg000940l", "utg000057l", "utg000941l", "utg000166l", "utg000782l", "utg000239l", "utg000833l", "utg000004l", "utg000856l", "utg000099l", "utg000719l", "utg000224l", "utg001056l", "utg000140l", "utg001028l", "utg000027l", "utg000868l", "utg000126l", "utg000822l", "utg000152l", "utg001005l", "utg000136l", "utg000763l", "utg000178l", "utg000806l", "utg000403l", "utg000759l", "utg000458l", "utg000757l", "utg000466l", "utg000753l", "utg000438l", "utg000971l", "utg000981l", "utg000002l", "utg000290l", "utg000818l", "utg000108l", "utg000667l"};
  //	vector<string> hap2 = {"utg000840l", "utg000503l", "utg001040l", "utg000423l", "utg000761l", "utg000546l", "utg000786l", "utg000925l", "utg000859l", "utg000499l", "utg000880l", "utg000596l", "utg001042l", "utg000372l", "utg000688l", "utg000655l", "utg000672l", "utg000666l", "utg001053l", "utg000478l", "utg000896l", "utg000405l", "utg000962l", "utg000333l", "utg000931l", "utg000721l", "utg000011l", "utg000534l", "utg000750l", "utg000130l", "utg000723l", "utg000270l", "utg000862l", "utg000271l", "utg000956l", "utg000102l", "utg000965l", "utg000218l", "utg000919l", "utg000273l", "utg000905l", "utg000223l", "utg000756l", "utg000116l", "utg001010l", "utg000951l", "utg000929l", "utg000298l", "utg000746l", "utg000231l", "utg000748l", "utg000681l", "utg000214l", "utg000831l", "utg000660l", "utg000964l", "utg000687l", "utg000885l", "utg000631l", "utg000476l", "utg000853l", "utg000216l", "utg000911l", "utg000256l", "utg001011l", "utg000072l", "utg001020l", "utg000117l", "utg000629l", "utg000630l", "utg000828l", "utg000129l", "utg000870l", "utg000769l", "utg000997l", "utg000041l", "utg000979l", "utg000251l", "utg000758l", "utg000051l", "utg001038l", "utg000277l", "utg000895l", "utg000210l", "utg000846l", "utg000125l", "utg000889l", "utg000094l", "utg000933l", "utg000133l", "utg001034l", "utg000048l", "utg000804l", "utg000078l", "utg000734l", "utg000149l", "utg000926l", "utg000594l", "utg000285l", "utg000010l", "utg000969l", "utg000115l", "utg000879l", "utg000854l", "utg000321l", "utg000324l", "utg000330l", "utg000847l", "utg000327l", "utg000960l", "utg000796l", "utg000650l", "utg000452l", "utg000713l", "utg000522l", "utg000848l", "utg000344l", "utg000829l", "utg000398l", "utg000978l", "utg000357l", "utg000735l", "utg000349l", "utg000680l", "utg000656l", "utg000639l", "utg000572l", "utg000766l", "utg000564l", "utg000861l", "utg000422l", "utg000637l", "utg000543l", "utg000700l", "utg000448l", "utg000963l", "utg000545l", "utg000805l", "utg000450l", "utg000679l", "utg000531l", "utg000745l", "utg000878l", "utg000898l", "utg000729l", "utg000686l", "utg000508l", "utg000752l", "utg000319l", "utg001006l", "utg000510l", "utg000607l", "utg000424l", "utg000809l", "utg000396l", "utg001055l", "utg000587l", "utg000988l", "utg000380l", "utg000812l", "utg000973l", "utg000736l", "utg000383l", "utg000913l", "utg000864l", "utg000529l", "utg000473l", "utg001013l", "utg001021l", "utg000544l", "utg000996l", "utg000447l", "utg000906l", "utg000429l", "utg000790l", "utg000337l", "utg000815l", "utg001017l", "utg000390l", "utg000732l", "utg000024l", "utg001039l", "utg000162l", "utg000743l", "utg000220l", "utg000697l", "utg000161l", "utg000624l", "utg000185l", "utg000907l", "utg000019l", "utg000917l", "utg000080l", "utg000773l", "utg000074l", "utg000775l", "utg000614l", "utg000922l", "utg000869l", "utg000354l", "utg000404l", "utg000952l", "utg000502l", "utg000691l", "utg000459l", "utg000768l", "utg000501l", "utg000726l", "utg000395l", "utg000845l", "utg000519l", "utg000704l", "utg000389l", "utg000891l", "utg000413l", "utg000765l", "utg000985l", "utg000589l", "utg000453l", "utg000388l", "utg000783l", "utg000348l", "utg000839l", "utg000460l", "utg000820l", "utg000430l", "utg000808l", "utg000542l", "utg000943l", "utg000420l", "utg000657l", "utg000802l", "utg000305l", "utg000649l", "utg000375l", "utg000940l", "utg000416l", "utg000941l", "utg000432l", "utg000782l", "utg000580l", "utg000833l", "utg000377l", "utg000856l", "utg000408l", "utg000719l", "utg000355l", "utg001056l", "utg000314l", "utg001028l", "utg000332l", "utg000868l", "utg000433l", "utg000822l", "utg000384l", "utg001005l", "utg000489l", "utg000763l", "utg000481l", "utg000806l", "utg000037l", "utg000759l", "utg000458l", "utg000757l", "utg000267l", "utg000753l", "utg000438l", "utg000971l", "utg000981l", "utg000002l", "utg000565l", "utg000818l", "utg000108l", "utg000667l"};	
//	vector<int> hap1_switches = {13, 21, 22, 39, 74, 82, 97, 98, 113};
//	vector<int> hap2_switches = {12, 21, 22, 39, 48, 50, 74, 82, 113};
	
	//add the last position of haplotype to the block ends as well	
	hap1_switches.push_back(hap1.size()-1);
	Node curnode;
	Node nextnode;
	string sequence = "";
	vector<string> sequence_blocks;
	curnode = graph.getNode(my_raw_id(hap1[0]));
	sequence += curnode.node_seq;
	for (int i = 1; i < hap1.size(); i++ ){
		
		int nodeid = my_raw_id(hap1[i]);	
		nextnode = graph.getNode(nodeid);		
	//	curnode = graph.getNode(nodeid);
	//	nextnode = graph.getNode(my_raw_id(hap1[i+1]));
	//	sequence += curnode.node_seq;
		//find out direction to find the correct NodePos pair for offsets
		pair<NodePos, NodePos> to_from = graph.getEdge(curnode, nextnode);
		
		int offset = graph.offsets[to_from];
		sequence += nextnode.node_seq.substr(offset);
		cout << "i: " << i << " curnode: " << curnode.node_id << " nextnode: " << nextnode.node_id << " sequence length: " << sequence.size() << endl;		
		if (find(hap1_switches.begin(), hap1_switches.end(), i) != hap1_switches.end()) {
			//if it's the last block end			
			if (i == hap1_switches[hap1_switches.size()-1]) {
				sequence_blocks.push_back(sequence);
			}
			else {
				sequence_blocks.push_back(sequence);
				sequence = "";		
				curnode = graph.getNode(my_raw_id(hap1[i+1]));
				sequence += curnode.node_seq;
				i += 1;
			}
		}
	//	i++;
	}
	cout << "number of sequence blocks: " << sequence_blocks.size() << endl;
	int sum;
	for (auto& block: sequence_blocks) {
		cout << block.size() << ", ";	
		sum += block.size();
	}
	cout << endl;
	cout << "length sequence: "<< sum << endl;	
	
return(0);	
*/
	//outfile: take stem of graph file name
//	string infofile = gfafile.substr(0,gfafile.find(".gfa"));  
	string infofile = "test";  
	ofstream myfile;
	myfile.open (infofile+"-bubbleinfo.txt");
	for (auto& chain: graph.chains) {
		myfile << "chain id: " << chain.id << endl;
		for (auto& bubble: graph.getChain(chain.id).bubbles) {
			myfile << "bubble id: " << bubble.id << endl;
			myfile << "node id: " ;
			for (auto& node: bubble.getNodes())	{
				myfile << node.node_id << ",";
			}
			myfile << endl;
		}	
  	}
	myfile.close();	
	


	AlignmentReader alignmentreader;
	alignmentreader.readAlignmentfile(alignmentfile, graph);
	
	cout << "alignments read" << endl;	
	cout << "number of alignments " << alignmentreader.alignments.size() << endl;	

	int testchain = 6;	
	
	cout << "number of alignments["<<testchain <<"]: " << alignmentreader.alignments[testchain].size() << endl;
//	cout << "long alignments of chain "<<testchain << ": " << endl;	
//	for (int i = 0; i < alignmentreader.alignments[testchain].size(); i ++){
//		if (alignmentreader.alignments[testchain].at(i).size() >= 4) {
//		alignmentreader.alignments[testchain].at(i).print();
//		}
//	}

	unordered_map<int, vector<vector<int>>> pathToAlleles;
	//TODO: test node -> allele relationship afterwards
	
	//allele of node 1241:
//	cout << "allele of 1241 before " << graph.nodes[1241].allele << endl;
	pathToAlleles = ChainsToReadset(graph);	
//	cout << "allele of 1241 after " << graph.nodes[1241].allele << endl;	
	cout << "chain paths computed " << endl;
	cout << "no. of chain paths " << pathToAlleles.size() <<  endl;


	unordered_map<int, unordered_map<int,vector<vector<int>>>> chainpathToAlleles;
	chainpathToAlleles = ChainsToReadsetDetailed(graph);
	cout << "chain paths computed " << endl;
	cout << "no. of detailed chain paths " << chainpathToAlleles.size() <<  endl;
	
	
	cout << "chain 8: " << endl;
	graph.printChain(8);
	for (auto& path: pathToAlleles[8]) {
		cout << "paths of 8 " << endl;
		for (auto& el: path) {
			cout << el << ",";		
		}	
		cout << endl;
	}
	

	
	cout << "         chainpaths 6:         " << endl;
	for (auto& path: chainpathToAlleles[6]) {
			cout << "bubble" << path.first << endl;
			for (auto& el: path.second) {
				for (auto& node: el) {
					cout << node << ",";			
				}
			cout << endl;		
			}	
		cout << endl;
	}	
	
	/**
	compute edit distances
	**/	
	
/*	cout << "compute edit distance between bubbles: " << endl;
	string editdistfile = gfafile.substr(0,gfafile.find(".gfa"));  
	ofstream edfile;
	edfile.open(editdistfile+"-editdistance.txt");

	map<int, map<int, vector<float>>> chain_to_editdist;
	for (auto& chainmap: chainpathToAlleles) {
		int chainid = chainmap.first;
		edfile << "chain id: " << chainid << endl;
		if (chainid == 6) {
			map<int, vector<float>> bubble_to_editdist;
			for (auto& bubblemap: chainmap.second) {
				
				int bubbleid = bubblemap.first;
				edfile << "bubble id: " << bubbleid << endl;
			//	if (bubbleid == 30 || bubbleid == 31 || bubbleid == 32 || bubbleid == 33 || bubbleid == 34 || bubbleid == 35 || bubbleid == 36 || bubbleid == 37 ||bubbleid == 38 ) {
				if (bubbleid == 12 || bubbleid == 13) { // || bubbleid == 14 || bubbleid == 15 || bubbleid == 16) {
					cout << "bubbleid: " << bubbleid << endl;
				vector<int> allelepath;
				int allele;
				vector<string> alleles;
				for (int it =0; it < bubblemap.second.size(); it++) {
					//every path is an allele
					edfile << "allele : " << it << endl;
					allelepath = bubblemap.second.at(it);
					string allelesequence = "";
					for (auto& node: allelepath) {
						edfile << node << ",";
						allelesequence += graph.getNode(node).node_seq;
					}
					edfile << endl;
					allele = it;
					alleles.push_back(allelesequence);
				}
				vector<float> eds;
				vector<int> eds_absolute;
				for (int i = 0; i < alleles.size(); i++){
					for (int j = 1; j < alleles.size(); j++) {
						if (i != j) {
						int editdist = EditDistDP(alleles[i], alleles[j]);
						float edit_ratio = -1;
						if (alleles[i].size() >= alleles[j].size())
							edit_ratio = ((float)editdist/(float)alleles[i].size())*100;
						else
							edit_ratio = ((float)editdist/(float)alleles[j].size())*100;
					//	eds.push_back(editdist);
					cout << "editdist: " << editdist << endl;
					edfile << "editdist: " << editdist << endl;
					cout << "edit_ratio in percent: " << to_string(edit_ratio) <<  endl;
					edfile << "edit_ratio in percent: " << to_string(edit_ratio) <<  endl;					
					cout << "size of nodes: " << alleles[i].size() << ", " << alleles[j].size() << endl;
					edfile << "size of nodes: " << alleles[i].size() << ", " << alleles[j].size() << endl;
						eds_absolute.push_back(editdist);
						eds.push_back(edit_ratio);
						}					
					}				
				}

			}
		}
			chain_to_editdist[chainid] = bubble_to_editdist;	
		}
	}
	edfile.close();
	
	cout << "edit distance between bubbles computed" << endl;

	return(0);
*/

	cout << "before alignmentsToReadset: " << alignmentreader.alignments.size() << endl;
	string readsetfile = alignmentfile.substr(0,alignmentfile.find(".gaf"));
	alignmentsToReadset(alignmentreader, graph, chainpathToAlleles, readsetfile);	
	cout << "after alignmentsToReadset: " << alignmentreader.alignments.size() << endl;
	


/*
TODO: overload == operator for chains
	for (auto& c: graph.chains) {
		vector<Node> c_nodes = c.listChain();
		sort(c_nodes.begin(), c_nodes.end());
		int equal = 0;
		for (auto& g: graph.chains) {
			vector<Node> g_nodes = g.listChain();
			sort(g_nodes.begin(), g_nodes.end());
    		if (c_nodes == g_nodes) {
				equal += 1;
				cout << "equal: " << c.id << " and " << g.id << endl;    
						
    		}		
		}	
	}
*/	
	return 0;
}



