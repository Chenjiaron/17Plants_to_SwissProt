# open file
fn = "pfam_info.tsv"
f = open(f"{fn}")

# data
rep_pfams = {}

# read file
while True:
    line = f.readline().strip()
    
    if not line:
        break
        
    tokens = line.split()
    hit = int(tokens[0])
    repId = tokens[1]
    pfams = set(tokens[2].split(';')[:-1])
    
    if not rep_pfams.get(repId):
        rep_pfams[repId] = []
    
    rep_pfams[repId].append([pfams, hit])

# close file
f.close()

# count the diversity of pfam
def compute_sum_hits(pfam_hits):
    L = len(pfam_hits)
    
    sum_hits = 0
    for i in range(L):
        [pfam, hits] = pfam_hits[i]
        
        sum_hits += hits
        
    return sum_hits

import math

rep_cov = {}
iteration = math.inf

check_id = 'AF-A0A009G5I8-F1-model_v3.cif'

for repId, pfam_hits in rep_pfams.items():
    
    # if there is only one hit at a cluster, we don't measure the consistency value
    if compute_sum_hits(pfam_hits) < 2:
        continue
        
    rep_cov[repId] = 0
    N = len(pfam_hits)    
    
    repId_hits = 0
    
    for i in range(N):
        pairwise_score = 0
        
        query_pfams = pfam_hits[i][0]
        query_hits = pfam_hits[i][1]
        repId_hits += query_hits
        
        query_N_pfams = len(query_pfams)
        
        for j in range(N):
            coverage = 0
            target_pfams = pfam_hits[j][0]
            target_hits = pfam_hits[j][1]
            
            # w/o self-pair
            if i == j :
                target_hits -= 1
            
            for pfam in query_pfams:
                if pfam in target_pfams:
                    coverage += 1
            
            coverage = coverage/ query_N_pfams * target_hits
            pairwise_score += coverage
        
        pairwise_score *= query_hits
        rep_cov[repId] += pairwise_score
    
    rep_cov[repId] /= (repId_hits**2 - repId_hits)

####
fon = "pfam-consistency_clusterId_cov.tsv"
fo = open(f'{fon}', 'w')

for key, value in rep_cov.items():
    fo.write(f'{key}\t{value}\n')

fo.close()
