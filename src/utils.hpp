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
