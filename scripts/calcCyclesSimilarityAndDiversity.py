from __future__ import print_function
import sys

def getCycles(fn):
    cycles=dict()
    for line in open(fn):
        rxes=[ x for x in line.split(';') if '>>' in x]
        if not len(rxes) in cycles:
            cycles[ len(rxes) ] = []
        inCycleCompounds=[]
        for rxpos in  range( len(rxes) ):
            secondRx=rxes[rxpos]
            firstRx=rxes[rxpos-1]
            firstRxProd=firstRx.split('>>')[1].split('.')
            secondRxSbs=secondRx.split('>>')[0].split('.')
            inCycleReagent = [ x for x in firstRxProd if x in secondRxSbs]
            if len(inCycleReagent) != 1:
                print("something wrong happend!")
                raise
            inCycleCompounds.append( inCycleReagent[0] )
        cycles[ len(rxes) ].append( inCycleCompounds)
    #print( [ (x, len(cycles[x])) for x in cycles] )
    return cycles


def calcDiversityOfCycles(cycles):
    #cycles is list of list
    allCompounds=[]
    for cycle in cycles:
        allCompounds.extend(cycle)
    uniqCmd=set(allCompounds)
    return len(uniqCmd)/ float(len(allCompounds) ), len(allCompounds), len(uniqCmd)


def calcSimilarityOfCycles(cycles):
    allCompounds=[]
    for cycle in cycles:
        allCompounds.extend(cycle)
    uniqCmd=set(allCompounds)
    mostFreqCmd=''
    mostFreq=0
    for cmd in uniqCmd:
        numCmd= allCompounds.count(cmd)
        if numCmd > mostFreq:
            mostFreq=numCmd
            mostFreqCmd=cmd
    hasMostCommon=[ cycle for cycle in cycles if (mostFreqCmd in cycle) ]
    return len(hasMostCommon) / float( len(cycles) )


if __name__ == "__main__":
    cycles=getCycles(sys.argv[1])
    print('cycle length', 'number of cycles', 'all compounds', 'unique compounds', 'diversity', 'similarity', sep=';')
    for cycleLen in cycles:
        diver, lenAllCmd, lenUniqCmd = calcDiversityOfCycles( cycles[ cycleLen])
        simil = calcSimilarityOfCycles(cycles[ cycleLen])
        print( cycleLen, len(cycles[cycleLen]), lenAllCmd, lenUniqCmd, diver, simil, sep=';')