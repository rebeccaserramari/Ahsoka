#include <iostream>
#include <sstream>

#include "argumentparser.hpp"
#include "timer.hpp"

#include "graph.hpp"
#include "alignmentreader.hpp"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cassert>

#include "/home/rebecca/work/whatshap-code/whatshap/src/read.h"
#include "/home/rebecca/work/whatshap-code/whatshap/src/readset.h"
#include "/home/rebecca/work/whatshap-code/whatshap/src/polyphase/readscoring.h"
#include "/home/rebecca/work/whatshap-code/whatshap/src/polyphase/trianglesparsematrix.h"
#include "/home/rebecca/work/whatshap-code/whatshap/src/polyphase/clustereditingsolver.h"
#include "/home/rebecca/work/whatshap-code/whatshap/src/polyphase/clustereditingsolution.h"
#include "/home/rebecca/work/whatshap-code/whatshap/src/polyphase/haplothreader.h"

//#include "stdlogger.h"
#include "filelogger.h"
#include <mutex>
using namespace std;

void get_coverage(ReadSet& readset, ClusterEditingSolution& clustering, vector<uint32_t> pos_index, vector<map<double, double>>& cov);
void get_pos_to_clusters_map(vector<map<double, double>>& coverage, uint32_t ploidy, vector<vector<GlobalClusterId>>& covMap, map<int,int>);
void add_to_vector(vector<map<double,double>>& cov);
void get_local_cluster_consensus(ReadSet& readset, ClusterEditingSolution& clustering, vector<vector<GlobalClusterId>>& cov_map, vector<map<double,double>>& consensus );
void get_single_cluster_consensus_frac(ReadSet& readset, ClusterEditingSolution& clustering, map<double,vector<uint32_t>>& rel_pos, vector<map<double, pair<double,double>>>& clusterwise_consensus);
vector<unsigned int> sort_asc(map<double, double>& M);
bool is_subset(vector<int> allelepath, AlignmentPath alignment, bool take_partial);
int getIndex(vector<unsigned int>* vec, unsigned int val);

void print(vector<int> vec){
	for (auto& el: vec) {
	cout << el << ",";
	}
	cout << endl;
}

template <typename T>
bool vector_contains(vector<T>* vec, T elem){
    bool result = false;
    for (auto & x : *vec){
        if (x == elem){
            result = true;
            break;
        }
    }
    return result;
}

void alignmentsToReadset(AlignmentReader& alignmentreader, Graph& graph, unordered_map<int, unordered_map<int, vector<vector<int>>>>& pathToAlleles, string readsetfile, bool shell_logging, vector<pair<int, int>>& size_sorting,std::mutex& g_display_mutex) {
	//std::mutex g_display_mutex;
	std::lock_guard<std::mutex> guard(g_display_mutex);
	Logger logger;
	if (shell_logging) {
		 //logger = make_shared<StdLogger>();
		 logger = make_shared<FileLogger>();
	}
	else {
		logger = make_shared<FileLogger>();
	}
	unordered_map<int, ReadSet*> newreadsets;
	
	int testchain = 0;	

	ofstream full_output;
	string tmp = readsetfile +"-result.txt";
	full_output.open(tmp,std::ios_base::app);


	for (auto& size: size_sorting) {
		auto chainmap = make_pair(size.second, pathToAlleles[size.second]);
	//for (auto& chainmap: pathToAlleles) {
		ReadSet* partial_readset;
		partial_readset = new ReadSet();
		int chainid = chainmap.first;
		logger->log_info("chain id: " + to_string(chainid));
		full_output << "chain id: " << chainid << endl;
		full_output << "size of chain: "<< chainmap.second.size() << endl; 

	//	if (chainid == testchain) {
		if (chainmap.second.size() > 1) {
		Read* read;
		ReadSet* readset;
		readset = new ReadSet();
		for (auto& bubblemap: chainmap.second) {
			int bubbleid = bubblemap.first;
			//logger->log_info("bubbleid, readset size: " + to_string(bubbleid) +" " + to_string(readset->size()));
			//logger->log_info("size bubblemap second " + to_string(bubblemap.second.size()));
			read = new Read(to_string(bubbleid),30,0,0,-1,"");
			vector<int> allelepath;
			int allele;
			map<string,string> readid_to_name;
			for (int it =0; it < bubblemap.second.size(); it++) {
				//every path is an allele
				allelepath = bubblemap.second.at(it);
				//print(allelepath);
				allele = it;
				//TODO: map here position(bubbleid) to allele to allelepath for the result?
				logger->log_info("number of alignments for chain " +to_string(chainid)+": "+to_string(alignmentreader.alignments[chainid].size()));
				bool read_exists = false;
				for (int al_it=0; al_it < alignmentreader.alignments[chainid].size(); al_it ++) {
					AlignmentPath alignment = alignmentreader.alignments[chainid].at(al_it);
					
					//if allelepath is a subset of alignment
					//alignment is AlignmentPath
					if (is_subset(allelepath, alignment, false)) {
				//	if (is_subset(allelepath, alignment, true)) {
						read_exists = true;
						string readid = alignment.name;
						//if the read does not exist yet in the readset, add it
						if (readset->getByName(readid, 0) == 0) {
							read = new Read(readid, alignment.id *100,0,0,-1,"");
							read->addVariant(bubbleid ,allele, 30);
							read->sortVariants();
							readset->add(read);
						}
						//if it already exists, add the variant
						else {
							//don't add variant if it already exists in read
							unordered_set<unsigned int> readpositions;
							readset->getByName(readid,0)->addPositionsToSet(&readpositions);
							if (readpositions.find(bubbleid) == readpositions.end()) {
								readset->getByName(readid,0)->addVariant(bubbleid ,allele, 30);
								readset->getByName(readid,0)->sortVariants();
							}						
						}					
					}				
				}
			}		
		}
		newreadsets[chainid] = readset;	
//		break; //for testing
//	}
	
	
	logger->log_info("               new readset               ");
	logger->log_info(" number of new readsets: " +to_string(newreadsets.size()));
	logger->log_info("number of readsets for chain "+to_string(chainid)+": "+to_string(newreadsets[chainid]->size()));	
	logger->log_info("print readsets for chain "+to_string(chainid)+": "+newreadsets[chainid]->toString());	
	
	ReadSet full_testset;
	full_testset = *newreadsets[chainid];

//	Read* read;
	IndexSet indices;
	for (int i = 0; i < full_testset.size(); i++) {
		read = full_testset.get(i);

		vector<int> mapqs;
		mapqs = read->getMapqs();
		assert(mapqs.size() == 1);
	//	if (read->getVariantCount() > 1 && mapqs.at(0) >= 97) {
		if (read->getVariantCount() > 1 && mapqs.at(0) >= 93) {
			indices.add(i);		
		}	
	}
	logger->log_info("indices computed");
	ReadSet* testset = full_testset.subset(&indices);	
	logger->log_info("testset size: "+to_string(testset->size()));
	logger->log_info("testset: " +testset->toString());


	/**
		Adding partial alignments for positions that are not fully covered
	
	**/

	//find boundary positions
	unordered_set<int> first_readpositions;
	unordered_set<int> last_readpositions;

	for (int i = 0; i < testset->size(); i++) {
		read = testset->get(i);
		first_readpositions.insert(read->firstPosition());
		last_readpositions.insert(read->lastPosition());

	}
	set<int> boundaries;
	for (auto& el: last_readpositions) {
		if (find(first_readpositions.begin(),first_readpositions.end(), el) == first_readpositions.end()) {
			boundaries.insert(el);		
			boundaries.insert(el+1);		
		}	
	}
	
	//add also gaps to the boundaries
	set<int> gaps;
	for (unsigned int i = 0; i < (*full_testset.get_positions()).back(); i++) {
		//if index is not in readset positions, it's a gap -> add to boundaries
		if (!vector_contains(testset->get_positions(), i)) {			
			gaps.insert(i);		
		}	
	}
	

	set<int> to_be_added;
	to_be_added.insert(boundaries.begin(), boundaries.end());
	to_be_added.insert(gaps.begin(), gaps.end());
	

	for (unsigned int i = 0; i < (*full_testset.get_positions()).back(); i++) {
		to_be_added.insert(i);					
	}
	logger->log_info("number of positions to be added: " + to_string(to_be_added.size()));
	for (auto& boundary: to_be_added) {

			vector<int> allelepath;
			int allele;
			int bubbleid = boundary;
	
		for (int it =0; it < pathToAlleles[chainid][boundary].size(); it++) {
				//every path is an allele
				allelepath = pathToAlleles[chainid][boundary].at(it);
				allele = it;
				
				for (int al_it=0; al_it < alignmentreader.alignments[chainid].size(); al_it ++) {
					AlignmentPath alignment = alignmentreader.alignments[chainid].at(al_it);
					
					//if allelepath is a subset of alignment
					//alignment is AlignmentPath
					if (is_subset(allelepath, alignment, true)) {
						string readid = alignment.name;
						//leave out alignments which start only in the overlap of the node
						//TODO: do not leave out the complete read, only the allele in the beginning! Idea: move this to the is_subset part: if allelepath is at beginning of alignment...
						//TODO: add the same for the last allele of the alignment						
						
						//if the read does not exist yet in the readset, add it
						if (partial_readset->getByName(readid, 0) == 0) {
							read = new Read(readid, alignment.id *100,0,0,-1,"");
							read->addVariant(bubbleid ,allele, 30);
							read->sortVariants();
							partial_readset->add(read);
							}
						//if it already exists, add the variant
						else {
							
							//don't add variant if it already exists in read
							unordered_set<unsigned int> readpositions;
							partial_readset->getByName(readid,0)->addPositionsToSet(&readpositions);
							if (readpositions.find(bubbleid) == readpositions.end() && (alignment.id*100) > 90) {
								
								partial_readset->getByName(readid,0)->addVariant(bubbleid ,allele, 30);
								partial_readset->getByName(readid,0)->sortVariants();
							}						
						}					
					}				
				}
		}	
	}

/**
	also subset the partial readset
**/
	logger->log_info("partial readset size: " +to_string(partial_readset->size()));
	logger->log_info("print partial readset before for chain "+to_string(chainid)+": " +partial_readset->toString());	

	Read* partial_read;
	IndexSet partial_indices;
	for (int i = 0; i < partial_readset->size(); i++) {
		partial_read = partial_readset->get(i);
		vector<int> mapqs;
		mapqs = partial_read->getMapqs();
		assert(mapqs.size() == 1);
	//	if (read->getVariantCount() > 1 && mapqs.at(0) >= 97) {
		if (partial_read->getVariantCount() > 1 && mapqs.at(0) >= 93) {
			partial_indices.add(i);		
		}	
	}
	ReadSet* partial_testset = partial_readset->subset(&partial_indices);	
/**
	partial readset subsetting ends
**/

	if (partial_testset->size() == 0) {
		logger->log_warning("No reads in ReadSet for chain " + to_string(chainid)+"!");
		continue;	
	}
	
	ofstream rsfile;
	string out_readsetfile = readsetfile +"-chain"+to_string(chainid)+"-readset.txt";
	rsfile.open(out_readsetfile);
	rsfile << "readsets for chain "<< chainid << ": " << newreadsets[chainid]->size() << endl;	
	rsfile << "print readsets for chain "<< chainid << ": " << newreadsets[chainid]->toString() << endl;	
	rsfile << "testset size: " << testset->size() << endl;
	rsfile << "testset: " << testset->toString() <<endl;
	rsfile << "partial testset size: " << partial_testset->size() << endl;
	rsfile << "partial testset: " << partial_testset->toString() << endl;
	rsfile.close();		
	
	//add the new partial testset to the testset
	testset = partial_testset;
	testset->sort();	
	
	ofstream finalrsfile;
	string finalreadsetfile = readsetfile + "-chain"+to_string(chainid)+"-readset_final.txt";	
	finalrsfile.open(finalreadsetfile);
	finalrsfile << "readset size: " << testset->size() << endl;
	finalrsfile << testset->toString() << endl;
	finalrsfile.close();
	
	int ploidy = 2;
	
	ReadScoring readscoring;
	TriangleSparseMatrix sim;
	//parameters: min_overlap and ploidy
	readscoring.scoreReadsetLocal(&sim, testset, 1,ploidy);
	ClusterEditingSolver clustereditingsolver(sim, false);
	ClusterEditingSolution clustering;
	clustering = clustereditingsolver.run();
	logger->log_info("number of clusters: "+to_string(clustering.getNumClusters()));
	vector<unsigned int>* pos;
	pos = testset->get_positions();

//	HaploThreader haplothreader(ploidy,32.0,8.0,true, 0);
	HaploThreader haplothreader(ploidy,32.0,8.0,false, 0);

	vector<Position> blockStarts; 
	vector<vector<GlobalClusterId>> covMap;
	vector<vector<double>> coverage; 
	vector<vector<uint32_t>> consensus;
	vector<unordered_map<uint32_t, uint32_t>> genotypes;
	vector<vector<GlobalClusterId>> path; 
	Position start = 0;
	Position end = pos->size();
	blockStarts.push_back(start);
	blockStarts.push_back(end);

	int numclusters = clustering.getNumClusters();
	int numpos = pos->size();
	vector<uint32_t> pos_index;
	for (auto& posi : *pos) {
		pos_index.push_back(posi);
	}
	
	numpos += 1;
	for (int i=0; i<numpos; i++){
		unordered_map<uint32_t, uint32_t> geno ( {{0,1},{1,1}} );
		genotypes.push_back(geno);	
	}
	
	vector<unsigned int>* readset_positions = testset->get_positions();	
	
	int lastpos = readset_positions->back();

	//map indices to those positions which are actually covered by reads
	map<int,int> indexmap;
	for (int i =0; i < lastpos+1; i++){
		if (find(readset_positions->begin(), readset_positions->end(), i) != readset_positions->end()) 
			indexmap[i] = i;
		else{
			indexmap[i] = -1;		
		}	
	}	
		
	logger->log_info("computing coverage: ");
	vector<map<double, double>> new_coverage(lastpos+1);

	get_coverage(*testset, clustering,pos_index, new_coverage );

	logger->log_info("computing coverage map: ");
	vector<vector<GlobalClusterId>> coverage_map;
	get_pos_to_clusters_map(new_coverage, ploidy, coverage_map, indexmap);

	logger->log_info("computing consensus: ");
	vector<map<double,double>> new_consensus;
	get_local_cluster_consensus(*testset, clustering, coverage_map, new_consensus );
	
	for (int i =0; i< new_consensus.size(); i++) {
		map<double,double> el = new_consensus.at(i);
	}
	

	//transform new_coverage from vector<map<double,double>> to vector<vector<double>>
	vector<vector<double>> new_coverage_vector(testset->get_positions()->size(), vector<double>());
	for (int i = 0; i < new_coverage.size(); i++) {
		if (new_coverage.at(i).size() == 0) {
			continue;
		}
		else {
			for (auto& map_entry: new_coverage.at(i)) {
				double c_id = map_entry.first;
				//new_coverage_vector.at(i)[c_id] = new_coverage.at(i)[c_id];
				int index = getIndex(testset->get_positions(), i);
				new_coverage_vector.at(index).push_back(new_coverage.at(i)[c_id]);
			}	
		}
	}	

	//transform new_consensus from vector<map> to vector<vector>	
	assert(new_consensus.size() == testset->get_positions()->size());
	vector<vector<uint32_t>> new_consensus_vector(testset->get_positions()->size(), vector<uint32_t>());
	for (int i = 0; i < new_consensus.size(); i++) {
		for (auto& map_entry: new_consensus.at(i)) {
			double c_id = map_entry.first;
			new_consensus_vector.at(i).push_back(new_consensus.at(i)[c_id]);
		}	
	}	

	assert(coverage_map.size() == new_coverage_vector.size() && coverage_map.size() == new_consensus_vector.size());

	//computing paths through the clusters	
	logger->log_info("computing paths: ");
	path = haplothreader.computePaths(start, end, coverage_map, new_coverage_vector, new_consensus_vector,genotypes);
	logger->log_info("path length: " +to_string(path.size()));	
	
	ofstream resfile;
	string resultfile = readsetfile + "-chain"+to_string(chainid)+"-result.txt";
	resfile.open(resultfile);
	
	vector<vector<uint32_t>> haps;
	for (int i=0; i< ploidy; i++) {
		vector<uint32_t> hap;
		vector<vector<int>> haplo;
		for (int j=0; j<path.size(); j++ ) {
			double c_id = path[j][i];
		//	uint32_t cons = new_consensus_vector[j][c_id];
			uint32_t cons = new_consensus[j][c_id];
			hap.push_back(cons);		
			vector<int> allelepath = pathToAlleles[chainid][testset->get_positions()->at(j)].at(cons);
			haplo.push_back(allelepath);
		}
		haps.push_back(hap);	
		set<int> usednodes;
		string logger_string = "";
		full_output << "haplotype " << i << ":" << endl;
		for (int len=0; len<haplo.size();len++) {
			auto node = haplo[len];
			
			for (int ind=0;ind<node.size()-1;ind++) {
				auto single = node[ind];
				auto next = node[ind+1];
				if (usednodes.count(single) != 0) {
					continue;				
				}			
				else {
					Node first = graph.getNode(single);
					Node sec = graph.getNode(next);
					pair<DirectedNode,DirectedNode> tup = graph.getEdge(first,sec);
					string firstdir;
					string secdir;
					(tup.first.end == 1)? firstdir = '+' : firstdir = '-';
					(tup.second.end == 1)? secdir = '+' : secdir = '-';
					logger_string += to_string(single)+'('+firstdir+')'+",";
					resfile << single	<< '(' << firstdir << ')' << ",";
					full_output << single	<< '(' << firstdir << ')' << ",";
					usednodes.insert(single);			
				}
			}
		}
		logger->log_info(logger_string);
		resfile << endl;
		full_output << endl;
		logger_string = "";
		for (int len=0; len<haplo.size();len++) {
			auto node = haplo[len];
			for (int ind=0;ind<node.size()-1;ind++) {
				auto single = node[ind];
				auto next = node[ind+1];
				Node first = graph.getNode(single);
				Node sec = graph.getNode(next);
				pair<DirectedNode,DirectedNode> tup = graph.getEdge(first,sec);
				string firstdir;
				string secdir;
				(tup.first.end == 1)? firstdir = '+' : firstdir = '-';
				(tup.second.end == 1)? secdir = '+' : secdir = '-';
				logger_string = to_string(tup.first.id) +','+firstdir +'|' +to_string(tup.second.id) + ',' + secdir + " ; ";			
			}
				
		}
		logger->log_info(logger_string);	
	}	
	resfile.close();
	
	for (auto& el: haps){
		cout << "hap: " << endl;
		for (int i = 0; i < el.size(); i++) {
			auto node = el[i];	
			cout << node <<"(" << pos->at(i) << ")" <<",";		
		}
		cout << endl;		
	}
	
	}
	
}
full_output.close();
return;
}

bool is_subset(vector<int> allelepath, AlignmentPath alignment, bool take_partial) {
	vector<int> alignment_raw = alignment.getRawIds();
	sort(alignment_raw.begin(), alignment_raw.end());


	
//	return(includes(alignment_raw.begin(), alignment_raw.end(), allelepath.begin(), allelepath.end()));	
	
	//return(middle_contained || full_contained);	
	if (!take_partial) {
		sort(allelepath.begin(), allelepath.end());
		bool full_contained = includes(alignment_raw.begin(), alignment_raw.end(), allelepath.begin(), allelepath.end()); 
		return(full_contained);	
	}
	else {
		allelepath.pop_back();
		allelepath.erase(allelepath.begin(), allelepath.begin()+1);
		sort(allelepath.begin(), allelepath.end());
		bool middle_contained = false;
		
	//	bool middle_contained = includes(alignment_raw.begin(), alignment_raw.end(), allelepath.begin()+1, allelepath.end()-1);
	//TODO for complex bubbles, it is sufficient to contain part of the allele instead of the full path (?)
		/*bool middle_contained = false;
		if (allelepath.size() > 1) {
			
			for (auto& it: allelepath) {
				if (find(alignment_raw.begin(), alignment_raw.end(), it) != alignment_raw.end())	{
					middle_contained = true;				
				}
			}

		}
		else {
		middle_contained = includes(alignment_raw.begin(), alignment_raw.end(), allelepath.begin(), allelepath.end());
	}
	*/
	//test if allelepath is at the beginning of the alignment - in that case, take overlap into account
/*	if (allelepath.at(0) == alignment_raw.at(0)) {
		//if the unique sequence of the node ends before the alignment starts, the alignment starts in the overlap to the second node -> do not add		
		if (alignment.firstnodelength < alignment.startpos) {
			middle_contained = false;
		}
		else {
			middle_contained = includes(alignment_raw.begin(), alignment_raw.end(), allelepath.begin(), allelepath.end());		
		}
	}
	else {
		middle_contained = includes(alignment_raw.begin(), alignment_raw.end(), allelepath.begin(), allelepath.end());
	}
*/
	middle_contained = includes(alignment_raw.begin(), alignment_raw.end(), allelepath.begin(), allelepath.end());
		return(middle_contained);
	}	
}

void get_local_cluster_consensus(ReadSet& readset, ClusterEditingSolution& clustering, vector<vector<GlobalClusterId>>& cov_map, vector<map<double,double>>& consensus ) {
	//Returns a list which for every position contains a dictionary, mapping a cluster id to
   // its consensus on this position
	map<double,vector<uint32_t>> relevant_pos;
	int num_vars = readset.get_positions()->size();

	for (int pos = 0; pos < num_vars; pos++) {
		for (auto& c: cov_map[pos]) {
			relevant_pos[c].push_back(pos);		
		}	
	}
	
	vector<map<double, pair<double,double>>> clusterwise_consensus;
	get_single_cluster_consensus_frac(readset, clustering, relevant_pos, clusterwise_consensus);
	vector<map<double,pair<double,double>>> whole_consensus;
	for (int pos = 0; pos < num_vars; pos++) {
		map<double,pair<double,double>> newdict;
		for (auto& c: cov_map.at(pos)) {
			newdict[c] = clusterwise_consensus[c][pos];
		}
		whole_consensus.push_back(newdict);	
	} 
	for (auto& pos: whole_consensus) {
		map<double,double> result_map;
		for (auto& [c_id, value] : pos) {
			result_map[c_id] = value.first;	
		}
		consensus.push_back(result_map);	
	}	 
	 

	
}

int getIndex(vector<unsigned int>* vec, unsigned int val) {
	int index = -1;
	auto it = find(vec->begin(), vec->end(), val);    
    if (it != vec->end()){
       index = it - vec->begin();
     }
	return(index);
}

void get_single_cluster_consensus_frac(ReadSet& readset, ClusterEditingSolution& clustering, map<double,vector<uint32_t>>& rel_pos, vector<map<double, pair<double,double>>>& clusterwise_consensus) {
	for (int i=0; i < clustering.getNumClusters(); i++ ){
		map<double, pair<double,double>> cluster_consensus;
		const vector<StaticSparseGraph::NodeId>& cluster = clustering.getCluster(i);
		vector<uint32_t> relevant_pos = rel_pos[i];	
		
		//Count zeroes and one for every position
		
    	int num_var = readset.get_positions()->back()+1;
    	vector<map<double, double>> poswise_allelecount(num_var);
    	
    	for (auto& readid: cluster)	 {
			Read* read = readset.get(readid);
			unordered_set<unsigned int> positions;
			read->addPositionsToSet(&positions);
			for (auto& pos: positions) {
				int allele = read->getAlleleFromPos(pos);
				//convert read index to the corresponding index 
				int index = -1;
				index = getIndex(readset.get_positions(), pos);
				assert(index != -1);
				poswise_allelecount.at(index)[allele] += 1;
			}
		}
		
		//Determine majority allele
		for (auto& pos: relevant_pos) {
			//convert pos into readset position?
			if (poswise_allelecount.at(pos).size() > 0) {
				int max_allele =0;
				double max_count =0;
				double sum_count =0;
				for (auto& allele: sort_asc(poswise_allelecount.at(pos))) {
					double cur_count = poswise_allelecount[pos][allele];
					sum_count += cur_count;
					if (cur_count > max_count) {
						max_allele = allele;
						max_count = cur_count;					
					}				
				}
				cluster_consensus[pos] = make_pair(max_allele, max_count/sum_count);			
			}	
			else {
				cluster_consensus[pos] = make_pair(0, 1.0);			
			}	
		}
 
		clusterwise_consensus.push_back(cluster_consensus);

	}
}




void get_coverage(ReadSet& readset, ClusterEditingSolution& clustering, vector<uint32_t> pos_index, vector<map<double, double>>& coverage){
	size_t num_vars = pos_index.size();
	size_t num_clusters = clustering.getNumClusters();

	//get all positions of readset and for each position, create an empty map and add to coverage vector
	vector<unsigned int>* readset_positions = readset.get_positions();

	vector<double> coverage_sum(readset_positions->back()+1,0);

	for (auto& p: *readset_positions){	
		map<double,double> covmap;
		coverage.at(p) = covmap;
	}
	for (unsigned int i = 0; i < num_clusters; i++){
		const vector<StaticSparseGraph::NodeId>& cluster = clustering.getCluster(i);
		for (auto& readid: cluster)	 {
			Read* read = readset.get(readid);
			unordered_set<unsigned int> positions;
			read->addPositionsToSet(&positions);

			for (auto& pos: positions) {
				map<double,double> testmap = coverage.at(pos);
				if (coverage.at(pos).size() == 0) {
					coverage.at(pos)[i] = 0;				
				}
				coverage.at(pos)[i] += 1;
				coverage_sum.at(pos) += 1;			
			}				
		}
	}
	
	for (int i = 0; i < coverage.size(); i++) {
		for (auto& map_entry: coverage.at(i)) {
			double c_id = map_entry.first;
			coverage.at(i)[c_id] = coverage.at(i)[c_id]/coverage_sum.at(i);	
		}	
	}
}

bool cmp(pair<double, double>& a, pair<double, double>& b) {
    return a.second > b.second;
}

bool cmp_asc(pair<double, double>& a, pair<double, double>& b) {
    return a.second < b.second;
}

vector<unsigned int> sort(map<double, double>& M) { 
  
    // Declare vector of pairs 
    vector<pair<double,double> > A; 
    vector<unsigned int> result; 
  
    // Copy key-value pair from Map 
    // to vector of pairs 
    for (auto& it : M) { 
        A.push_back(it); 
    } 
  
    // Sort using comparator function 
    sort(A.begin(), A.end(), cmp); 
  
    // Print the sorted value 
    for (auto& it : A) { 
  			result.push_back(int(it.first));
    }
    return(result); 
} 

vector<unsigned int> sort_asc(map<double, double>& M) { 
  
    // Declare vector of pairs 
    vector<pair<double,double> > A; 
    vector<unsigned int> result; 
  
    // Copy key-value pair from Map 
    // to vector of pairs 
    for (auto& it : M) { 
        A.push_back(it); 
    } 
  
    // Sort using comparator function 
    sort(A.begin(), A.end(), cmp_asc); 
  
    // Print the sorted value 
    for (auto& it : A) { 
  			result.push_back(int(it.first));
    }
    return(result); 
}

void get_pos_to_clusters_map(vector<map<double, double>>& coverage, uint32_t ploidy, vector<vector<GlobalClusterId>> &covMap, map<int,int> indexmap) {

	/* 	For every position, computes a list of relevant clusters for the threading
	algorithm. Relevant means that the relative coverage is at least 1/8 of
	what a single haplotype is expected to have for the given ploidy. Apart
	from that, at least <ploidy> and at most <2*ploidy> many clusters are
	selected to avoid exponential blow-up.
	*/
	vector<GlobalClusterId> sorted_cids;
	for (int pos = 0; pos < coverage.size(); pos++) {
		if (indexmap[pos] == -1) {
			continue;		
		}
		else {
		    sorted_cids = sort(coverage[pos]);
		    size_t cut_off = min(uint32_t(sorted_cids.size()), 2*ploidy);
		    for (int i = ploidy; i < min(uint32_t(sorted_cids.size()), 2*ploidy);i++) {
				if (coverage[pos][sorted_cids[i]] < (1.0 / (8.0 * ploidy))) {
					cut_off = i;
					break;			
				}	    
		    }
		    vector<GlobalClusterId> cut_sorted_cids(cut_off);
			copy(sorted_cids.begin(), sorted_cids.begin()+cut_off, cut_sorted_cids.begin());
			covMap.push_back(cut_sorted_cids);
	 	}
	} 
	return;
}
