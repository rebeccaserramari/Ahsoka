import sys
import seaborn as sns
from matplotlib import pyplot as plt

bubblefile = sys.argv[1]
chain_to_bubblecount = {}
with open(bubblefile) as bfile:
	for i,line in enumerate(bfile):
		if line.strip().split(':')[0] == "chain id":
			chain_id = int(line.strip().split(':')[1])
		if line.strip().split(':')[0] == "bubble id":
			nextline = next(bfile)
			if chain_id not in chain_to_bubblecount:
				chain_to_bubblecount[chain_id] = 0 		
			chain_to_bubblecount[chain_id] += 1
			
print("number of chains: ", len(chain_to_bubblecount.keys()))			
print("maximum length: ", max(chain_to_bubblecount.values()))			
print("minimum length: ", min(chain_to_bubblecount.values()))			
print("number of singleton chains: ", len([i for i in chain_to_bubblecount.values() if i == 1]))
print("lengths of the longest 30 chains: ", sorted(chain_to_bubblecount.values(), reverse=True)[:30])

cutoff = 2
sns.histplot(x=[i for i in chain_to_bubblecount.values() if i >= cutoff])
plt.xlabel("bubble chain length")
plt.title("Distribution of lengths of bubble chains longer than " +str(cutoff-1))
plt.savefig(bubblefile.split('.txt')[0] + '-bubble_hist.pdf')
plt.show()
