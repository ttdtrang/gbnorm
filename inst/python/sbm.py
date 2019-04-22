import graph_tool as gt
import graph_tool.inference
import sys
import numpy as np
from timeit import default_timer as timer

if (len(sys.argv) < 3):
    print('''Usage:
    %s <graphFileName> <outPrefix> [nruns=3]
    ''' % sys.argv[0])
    exit()

graphFileName = sys.argv[1] # '../data/corr-graphs/g50.graphml'
outPrefix = sys.argv[2] # '../data/corr-graphs/g50.graphml'
outFileName = outPrefix + '.tsv' # '../data/corr-graphs/g50.graphml'
outReport = outPrefix + '.entropy' # '../data/corr-graphs/g50.graphml'
NRUNS = 3
if (len(sys.argv) >= 4): NRUNS = int(sys.argv[3])
g = gt.load_graph(graphFileName)
bLabels = np.zeros((g.num_vertices(), NRUNS),dtype=np.int8)
bScore = np.zeros(NRUNS)
t = timer()
for i in range(NRUNS):
    state = gt.inference.minimize_blockmodel_dl(g,deg_corr=False)
    bLabels[:,i] = state.get_blocks().get_array()
    bScore[i] = state.entropy()
dt1 = timer() - t
print("Input graph: %s" % graphFileName)
print("Number of vertices: %s" % g.num_vertices())
print("Number of runs: %s" % NRUNS)
print("Total time (seconds)  : %.3f" % dt1)
print("Average time (seconds): %.3f" % (dt1/NRUNS))
print("Best model (0-based index): %s" % bScore.argmin())

np.savetxt(outReport, bScore)

vname = g.vertex_properties['name']
with open(outFileName, 'w') as ofile:
    for v in g.vertices():
        labels = '\t'.join(bLabels[int(v),:].astype('str'))
        ofile.write('%s\t%s\t%s\n' % (v,vname[v],labels))

