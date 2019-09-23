import graph_tool as gt
from graph_tool import topology as tp
from multiprocessing import Process, Queue
import timeit


# search for cycles in graph using graph_tool library
# graph_tool is not very popular library however to my best knowledge it is the most 
# efficient lib for graph in sense of memory usage and speed
#
# adjust MAXPROC and MAXCYCLESIZE


idx=0
def buildGraph(dane):
    def _getNum(sm):
        global idx
        if not(sm in smiToNum):
            smiToNum[sm]=idx
            idx+=1
        return smiToNum[sm]
    ###
    graph=gt.Graph()
    #    for smi in allRes:
    #        rxSmi='.'.join(allRes[smi]['parents'])+'>>'+smi+':::'+str(allRes[smi]['rxid'])
    #        for subst in allRes[smi]['parents']:
    #            G.add_edge( subst, smi, rxinfo=rxSmi)
    #        if 'otherParents' in allRes[smi]:
    #            for oparents in allRes[smi]['otherParents']:
    #                rxSmi='.'.join(oparents['p'])+'>>'+smi+':::'+str(oparents['id'])
    #                for subst in oparents['p']:
    #                    G.add_edge(subst, smi, rxinfo=rxSmi)
    smiToNum=dict()
    numToSmi=dict()
    rxMap=dict()
    for smi in dane['results']:
        rxSmi='.'.join(dane['results'][smi]['parents'])+'>>'+smi+':::'+str(dane['results'][smi]['rxid'])
        for sbs in dane['results'][smi]['parents']:
            graph.add_edge( _getNum(sbs), _getNum(smi) )
            rxMap[(_getNum(sbs), _getNum(smi))]=rxSmi
        if 'otherParents' in dane['results'][smi]:
            for parents in dane['results'][smi]['otherParents']:
                rxSmi = '.'.join(parents['p'])+'>>'+smi+':::'+str(parents['id'])
                for sbs in parents['p']:
                    graph.add_edge( _getNum(sbs), _getNum(smi) )
                    rxMap[ (_getNum(sbs), _getNum(smi) )]=rxSmi
    return {'graph':graph, 'map':rxMap}


def getRxSmi(cykl, mapa):
    allRx=[]
    for i in range(1, len(cykl)):
       allRx.append( graphInfo['map'][(cykl[i-1],cykl[i])] )
    allRx.append(graphInfo['map'][(cykl[-1],cykl[0])] )
    return allRx

#epox=Chem.MolFromSmarts('[CX4]([O]1)[CX4]1')
def findPaths(qin, qout):
    while True:
        rx = qin.get()
        if rx == 'STOP':
            qout.put("END")
            return 0
        numPath=0
        MAXCYCLESIZE=70
        for path in tp.all_paths(graphInfo['graph'], rx[1], rx[0], cutoff=maxSizeOfCycle):
            thisCyc=frozenset(path)
            qout.put( thisCyc)
            numPath+=1
        print("CDONE", rx, numPath)
    return 1


if __name__ == "__main__":
    import pickle as cp
    import sys
    dane=cp.load( open(sys.argv[1], 'rb'), encoding='latin1' )
    graphInfo=buildGraph(dane)
    addedRx=set()
    MAXPROC=12
    allProcs = []
    qin = Queue(maxsize=123000)
    qout = Queue()
    for i in range(MAXPROC):
        p=Process( target=findPaths, args=(qin, qout) )
        p.start()
        allProcs.append( p)
    allSubmited=0
    for  rx in graphInfo['map']:
        qin.put(rx)
        allSubmited+=1

    print( "ALL IN Q", allSubmited)
    allAdded=set()
    cycStat={}
    sumCycles=0
    itr=0
    while True:
        cycle = qout.get()
        itr+=1
        if itr % 2000 ==0:
            print("ITR",itr, timeit.default_timer() )
        if not cycle in allAdded:
            allAdded.add(cycle)
        else:
            continue
        lenc = len(cycle)
        if not lenc in cycStat:
            cycStat[lenc]=0
        cycStat[lenc]+=1
        sumCycles +=1
        if sumCycles %100 == 0:
            print("STAT", cycStat)
            #sys.stdout.flush()
        #print("SC", sumCycles)
