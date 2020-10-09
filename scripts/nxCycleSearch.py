import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem

class rxGraph:
    def __init__(self, allRes, allRx=dict()):
        #self.allRes=allRes
        self.graph=self.buildGraph(allRes)
        self.allRx=allRx
        self.allSides=set()
        self.allNodes=[x for x in self.graph.nodes() ]

    def buildGraph(self, allRes):
        G = nx.DiGraph()
        for smi in allRes:
            rxSmi='.'.join(allRes[smi]['parents'])+'>>'+smi+':::'+str(allRes[smi]['rxid'])
            for subst in allRes[smi]['parents']:
                G.add_edge( subst, smi, rxinfo=rxSmi)
                #if 'CC(C)=O' in rxSmi:
                #    print "ADDED1", subst, smi, rxSmi
            if 'otherParents' in allRes[smi]:
                for oparents in allRes[smi]['otherParents']:
                    rxSmi='.'.join(oparents['p'])+'>>'+smi+':::'+str(oparents['id'])
                    for subst in oparents['p']:
                        G.add_edge(subst, smi, rxinfo=rxSmi)
                        #if 'CC(C)=O' in rxSmi:
                        #    print "ADDED", subst, smi, rxSmi
        return G

    def cycleWithRx(self, rx):
        self.graph.edge

if __name__ == "__main__":
    import sys
    from pprint import pprint
    import cPickle as cp
    allData=cp.load( open(sys.argv[1]) )
    allSmi=allData['results'].keys()
    graphObj=rxGraph( allData['results'] )
