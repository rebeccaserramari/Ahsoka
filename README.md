# Ahsoka
Haplotype assembly of diploid and polyploid genomes based on assembly graphs and long reads 


### Usage:

  Ahsoka phase -g \<graphfile\> -a \<alignmentfile\>

### Input:  
  
  \<graphfile\>: an assembly graph in .gfa format
  
  \<alignmentfile\>: contains alignments of long reads (advised: ONT ultra-long) to the graph in .gaf format
  
### Output:
  
  outputs a file in the following format:
  
    chain id: xx
    size of chain: xx (= number of bubbles)
    haplotype 0:
    node_1(+),node_2(-), ... node_n(+)
    haplotype 1:
    ...
    
    haplotype i:
    ...
    
 where i = <ploidy>-1 and (+) and (-) denote whether the node occurs in forward or in reverse direction.
