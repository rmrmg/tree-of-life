#import networkx as nx
#import math

from rdkit import Chem
lifeSmiles= set(['NCC(=O)O', 'CC(N)C(=O)O', 'NC(CC(=O)O)C(=O)O', 'NC(CO)C(=O)O', 'O=C(O)C1CCCN1', 'CC(O)C(N)C(=O)O', 'Nc1ncnc2[nH]cnc12', 'Nc1nc2[nH]cnc2c(=O)[nH]1', 'O=c1cc[nH]c(=O)[nH]1',
'Nc1cc[nH]c(=O)n1', 'O=c1[nH]cnc2[nH]cnc12', 'Nc1nc(=O)[nH]c2[nH]cnc12', 'O=c1[nH]c(=O)c2nc[nH]c2[nH]1', 'O=c1[nH]c(=O)c2[nH]c(=O)[nH]c2[nH]1', 'O=CC(O)CO', 'O=CC(O)C(O)CO',
'OCC1OC(O)C(O)C1O', 'O=CC(O)C(O)C(O)C(O)CO', 'OCC1OC(O)(CO)C(O)C1O', 'O=C(CO)C(O)C(O)C(O)C(O)CO', 'Nc1ncnc2c1ncn2C1OC(CO)C(O)C1O', 'Nc1nc2c(ncn2C2OC(CO)C(O)C2O)c(=O)[nH]1',
'O=c1ccn(C2OC(CO)C(O)C2O)c(=O)[nH]1', 'Nc1ccn(C2OC(CO)C(O)C2O)c(=O)n1', 'CC(O)C(=O)O', 'O=C(O)C=CC(=O)O', 'NC(N)=O', 'O=C(O)CC(O)C(=O)O', 'NC(CCC(=O)O)C(=O)O', 'Cc1c[nH]c(=O)[nH]c1=O',
'O=C(O)CCC(=O)O', 'OCC(O)CO', 'NC(Cc1ccccc1)C(=O)O', 'NCC(=O)NCC(=O)O', 'NCC(=O)NC(CC(=O)O)C(=O)O', 'NCC(=O)NC(CO)C(=O)O', 'CC(NC(=O)CN)C(=O)O', 'CC(N)C(=O)NCC(=O)O',
'CC(N)C(=O)NC(CC(=O)O)C(=O)O', 'CC(N)C(=O)NC(CO)C(=O)O', 'CC(N)C(=O)NC(C)C(=O)O', 'NC(CO)C(=O)NCC(=O)O', 'NC(CO)C(=O)NC(CC(=O)O)C(=O)O', 'NC(CO)C(=O)NC(CO)C(=O)O',
'CC(NC(=O)C(N)CO)C(=O)O', 'CSCCC(N)C(=O)O', 'NCC(=O)NCC(=O)NC(CC(=O)O)C(=O)O', 'NCC(=O)NCC(=O)NCC(=O)O', 'NCC(=O)NCC(=O)NC(CO)C(=O)O', 'CC(NC(=O)CNC(=O)CN)C(=O)O',
'NCC(=O)NCC(=O)NCC(=O)NCC(=O)O', 'NC(CO)C(=O)NCC(=O)NC(CC(=O)O)C(=O)O', 'NC(CO)C(=O)NCC(=O)NC(CO)C(=O)O', 'NC(CO)C(=O)NCC(=O)NCC(=O)O', 'CC(N)C(=O)NCC(=O)NCC(=O)NCC(=O)O',
'CC(NC(=O)CNC(=O)C(N)CO)C(=O)O', 'CC(N)C(=O)NCC(=O)NC(C)C(=O)O', 'CC(N)C(=O)NCC(=O)NC(CO)C(=O)O', 'CC(N)C(=O)NCC(=O)NC(CC(=O)O)C(=O)O', 'NC(CO)C(=O)NCC(=O)NCC(=O)NCC(=O)O',
'CC(N)C(=O)NCC(=O)NCC(=O)O', 'CC(NC(=O)C(N)Cc1ccccc1)C(=O)O', 'NC(Cc1ccccc1)C(=O)NCC(=O)O', 'NC(Cc1ccccc1)C(=O)NC(CC(=O)O)C(=O)O', 'NC(Cc1ccccc1)C(=O)NCC(=O)NCC(=O)O',
'NC(Cc1ccccc1)C(=O)NC(CO)C(=O)O', 'CC(NC(=O)CNC(=O)C(N)Cc1ccccc1)C(=O)O', 'CC(NC(=O)C(Cc1ccccc1)NC(=O)CN)C(=O)O', 'CC(N)C(=O)NCC(=O)NC(Cc1ccccc1)C(=O)O',
'NCC(=O)NC(Cc1ccccc1)C(=O)NCC(=O)O', 'CC(N)C(=O)NC(Cc1ccccc1)C(=O)NCC(=O)O', 'NC(CO)C(=O)NC(Cc1ccccc1)C(=O)O', 'NCC(=O)NC(Cc1ccccc1)C(=O)O', 'CC(N)C(=O)NC(Cc1ccccc1)C(=O)O',
'NCC(=O)NCC(=O)NC(Cc1ccccc1)C(=O)O', 'O=CC(=O)O', 'CC(=O)C(=O)O', 'O=C(O)C(=O)O', 'CC(=O)O', 'O=C(O)CC(=O)C(=O)O', 'O=C(O)CC(C(=O)O)C(O)C(=O)O', 'O=C(O)CC(O)(CC(=O)O)C(=O)O',
'O=C(O)C=C(CC(=O)O)C(=O)O', 'O=C(O)CC(=O)O', 'O=C(O)CCC(=O)C(=O)O'])


#lifeSmiles += lifeS


lifeSmiles=[ Chem.MolToSmiles(Chem.MolFromSmiles(s),True) for s in lifeSmiles]


condValue={'SA':-3, 'A':-2, 'WA':-1, 'N':0, 'WB':1, 'B':2, 'SB':3}
initialSubstrate=['C#N', 'N#N', 'C', 'N', 'O','S']


def getConValues(condList):
    newList =[]
    for cond in condList:
        cond= [ x for x in cond.split('+') if x != 'LA']
        
    return newList

def parseHistoryEntry(entry):
    condition=[]
    rxsplit=entry['smiles'].split('>>')
    rx={'substrates':rxsplit[0].split('.'), 'product':rxsplit[1]   }
    if 'tandem' in entry:
         condition=entry['tandem'][0]['condAB'] #.split('.')
    else:
        condition=entry['condAB'] #.split('.')
    try:
        #print "CONDITIO", condition
        condition = [condValue[x] for x in condition if not 'LA' in x]
        #condition= getConValues(condition )  #[ condValue[x] for x in condition if x != 'LA']
        if 0 in condition:
            if -2 in condition:
                if not -1 in condition:
                    condition.append(-1)
            if 2 in condition:
                if not 1 in condition:
                    condition.append(1)
        rx['condition']= tuple( sorted(condition) )
    except:
        print "\n PROBLERM WITH RX", entry
        #return
        #continue
        raise
    return rx

def makeCondVariant( dic):
    doUpTo=19999999
    iteration=0
    for smi in  dic['results']:
        ignoreit=False
        if smi in ['CSC(=O)C(C)OC(=O)CCN1CC(=O)NC1=O', 'CC(=O)C(C#N)NCC(=O)N(CCC(=O)O)CC(=O)O']:
            print "ignore PROBLEM"
            continue
        #if not( smi in lifeSmiles):
        #    continue
        #print "HIST:", dic['results'][smi]['history']
        
        allRx=[]
        for entry in dic['results'][smi]['history']:
            parsedEntry=parseHistoryEntry(entry)
            if parsedEntry:
                allRx.append( parsedEntry )
        #try:
        pathCond=getAllPath(allRx,smi)
        #except:
        #    print "\n\nFAILED ON", smi, "ALLRX", allRx,"\n"
        #    raise
        minPathLen=  [ pathCond['minLen'][0],  pathCond['minLen'][1]  , len(allRx), 0 ]
        minPathValue=[ pathCond['minValue'][0],  pathCond['minValue'][1], len(allRx), 0 ]
        alterMinLen=[]
        alterMinValue=[]
        #print "???",pathCond
        #raise
        minNumRx=len(allRx)
        minNumVariant=0
        #minRxNum=len(allRx)
        numOfRxPaths=1
        #print "....SMI ", smi, "numRx", len(allRx), "minLen", pathCond['minLen'], "minValue", pathCond['minValue'], "all posibible", len(pathCond['allPaths'])
        #print "ALLRX", allRx
        if 'alterHistories' in dic['results'][smi]:
            for alteridx, onealter in  enumerate(dic['results'][smi]['alterHistories']):
                print "smi", smi, "alternative", alteridx
                if alteridx> 1500:
                    ignoreit=True
                    break
                alteridx+=1
                allAlter=[]
                for entry in onealter:
                    parsedEntry = parseHistoryEntry(entry) 
                    if parsedEntry:
                        allAlter.append( parsedEntry)
                #print "SMI", smi, "ALTER", len(allAlter)
                pathCond=getAllPath(allAlter,smi)
                #print " .........", pathCond['minValue']
                numOfRxPaths+=1
                if minNumRx> len(allAlter):
                    minNumRx= len(allAlter)
                    minNumVariant=alteridx
                if pathCond['minValue'][0] < minPathValue[0]:
                    minPathValue=[ pathCond['minValue'][0], pathCond['minValue'][1], len(allAlter), alteridx ]
                    alterMinValue=[]
                elif pathCond['minValue'][0] == minPathValue[0]:
                    #print "PPP", pathCond
                    if len(pathCond['minValue'][1] ) < len(minPathValue[1]):
                        minPathValue=[ pathCond['minValue'][0], pathCond['minValue'][1], len(allAlter),alteridx ]
                    elif len(pathCond['minValue'][1] ) == len(minPathValue[1]):
                        if len(allAlter) < minPathValue[2]:
                            #print "UPDATE", pathCond, "UP MIN", minPathValue
                            #print "AA", len(allAlter) , "OLD", minPathValue[2]
                            minPathValue=[ pathCond['minValue'][0], pathCond['minValue'][1], len(allAlter), alteridx ]
                    alterMinValue.append( [pathCond['minValue'][0], pathCond['minValue'][1], len(allAlter),alteridx ] )

                if pathCond['minLen'][0] < minPathLen[0]:
                    minPathLen=[ pathCond['minLen'][0], pathCond['minLen'][1], len(allAlter), alteridx ]
                    alterMinLen=[]
                elif pathCond['minLen'][0] == minPathLen[0]:
                    oldPathValue=pathValue( minPathLen[1]  )
                    newPathValue=pathValue(pathCond['minLen'][1])
                    if newPathValue < oldPathValue:
                        minPathLen=[ pathCond['minLen'][0], pathCond['minLen'][1], len(allAlter), alteridx ]
                    elif newPathValue == oldPathValue:
                        if len(allAlter) < minPathValue[2]:
                            minPathLen=[ pathCond['minLen'][0], pathCond['minLen'][1], len(allAlter), alteridx]
                    alterMinLen.append( [ pathCond['minLen'][0], pathCond['minLen'][1], len(allAlter), alteridx] )
                #print "\t", "numRx", len(allAlter), "minLen", pathCond['minLen'], "minValue", pathCond['minValue'], "all posibible", len(pathCond['allPaths'])
        if ignoreit:
            print "IGNORE SMILES", smi
        life="NO"
        if smi in lifeSmiles:
            life="LIFE"
        #print life, "\t",  smi, "\t",minPathLen[0]-1, "\t", "MINLEN", minPathLen, "MINVALUE", minPathValue, "MIN rx len", minNumRx, "rx paths", numOfRxPaths
        #print life, '\t', smi, '\t', minPathLen[0]-1,  '\t', minPathLen[1][0], '\t', minPathLen[1][1], '\t', minPathLen[2], '\t',
        #print 'VALUE', '\t', minPathValue[0], '\t', minPathValue[1][0], '\t', minPathValue[1][1], '\t', minPathValue[2]
        try:
            condFromSubst=minPathLen[1][0][::-1]
            numOfRxFromSubst=minPathLen[1][1][::-1]
            fullPath=[]
            for jjj, repeat in enumerate(numOfRxFromSubst):
                for iii in range(repeat):
                    fullPath.append( condFromSubst[jjj])
            print life, minNumRx, '\t', smi, '\t', minPathLen[0]-1,  '\t', fullPath, '\t', minPathLen[2], minPathLen[3], '\t :::', alterMinLen, 
    
            condFromSubst=minPathValue[1][0][::-1]
            numOfRxFromSubst=minPathValue[1][1][::-1]
            fullPath=[]
            for jjj, repeat in enumerate(numOfRxFromSubst):
                for iii in range(repeat):
                    fullPath.append( condFromSubst[jjj])
            print '|||VALUE', '\t', minPathValue[0], '\t', fullPath, '\t', minPathValue[2], minPathValue[3], "\t :::", alterMinValue
            if life == 'LIFE':
                if minPathLen[3] == 0:
                    print "MIN PATH LEN", minPathLen[3], dic['results'][smi]['history']
                else: 
                    print "MIN PATH LEN", minPathLen[3], dic['results'][smi]['alterHistories'][minPathLen[3] -1]
                if minPathValue[3] == 0:
                    print "MIN PATH VALUE", minPathValue[3], dic['results'][smi]['history']
                else:
                    print "MIN PATH VALUE", minPathValue[3], dic['results'][smi]['alterHistories'][minPathValue[3] -1]
            if minNumRx != minPathLen[2] and minNumRx != minPathValue[2]:
                print "SHORTER more condiftion changes", minNumVariant
        except:
            print "PROBLEM WITH SMILES", smi
        if iteration > doUpTo:
            break
        iteration+=1

def getAllPath(allRxOryg,smi):
    ##
    readyPath=set()
    growingPath=[]
    ##get last step:
    lastStep=[rx for rx in allRxOryg if rx['product'] == smi]
    if len(lastStep) !=1:
        #print "\n no one last step", "ALL", allRxOryg
        return {'minLen':(99999, [-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,]), \
              'minValue':(99999, [-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,]), \
                "allPaths": set([ (-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,), ])}
        #raise
    lastStep=lastStep[0]
    for cond in lastStep['condition']:
        growingPath.append( { 'path':[], 'nextCondition':cond, 'rxToCheck':[lastStep,]} )
    #print "BEFORE WHILE", growingPath
    pfound=0
    while True:
        toGrow=[]
        for path in growingPath:
            #print "H", path
            allRetPaths=growPath( path, allRxOryg)
            #print "II", allRetPaths
            for (longerPath,status) in allRetPaths:
                if status == 'ready':
                    if longerPath:
                        #print "longerPath:", longerPath
                        readyPath.add( tuple( longerPath ) )
                    else:
                        #print "line 122 INP", smi,"  \nPATH:", path, "\nALL", allRxOryg
                        #raise
                        pass
                else:
                    toGrow.append( longerPath)
                    #print "\n---LP", longerPath
        if len(toGrow) == 0:
            break
        growingPath=toGrow[:]
        pfound+=1
        if pfound> 8400:
            print "TO MANY PATH for", smi
            break
        #for p in growingPath:
        #    #print "G", p
    minLen=[9999999, None]
    minValue=[9999999, None]
    alterMinLen=[]
    alterMinValue=[]
    #print "ALL", readyPath
    #raise
    for path in readyPath:
        pathVal=pathValue(path)
        ##find the shortest len (if len is equal chose this with lower changeValue if is equal chose this with lowest num of 
        if pathLen(path) < minLen[0]:
            minLen=[ pathLen(path), path]
            alterMinLen=[]
        elif pathLen(path) == minLen[0]:
            oldPathVal=pathValue( minLen[1])
            if oldPathVal> pathVal:
                minLen=[ pathLen(path), path]
            alterMinLen.append( [pathLen(path), path] )
            #elif oldPathVal == pathVal:
            #    rxNumOld=rxInPath(minLen[1] )
            #    rxNumNew=rxInPath(path)
            #    if rxNumNew < rxNumOld:
            #        minLen=[ pathLen(path), path]
            #        print "NEWER????"
            #    #elif rxNumNew == rxNumOld:
            #    #    print "EQUAL PATH", path, "MINLEN", minLen
        if pathVal < minValue[0]:
            minValue=[pathVal, path]
            alterMinValue=[]
        elif pathVal == minValue[0]:
            if pathLen(path) < pathLen( minValue[1] ):
                minValue=[pathVal, path]
            alterMinValue.append( [pathVal, path] )
    #print "   "
    #print "readyPATH", [len(p) for p in readyPath]
    #print "R", readyPath
    #print "MIN len", minLen, "min val", minValue
    if minLen[1] == None:
        #print "SMI", smi, "RP", readyPath, "\nALLRX", allRxOryg, "\n"
        #raise
        return {'minLen':(99999, [-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,]), \
              'minValue':(99999, [-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,]), \
                "allPaths": set([ (-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,-10,10,), ])}
        #return False
        #raise
    #print "XXXXXX", {'minLen':minLen, 'minValue': minValue, "allPaths": readyPath}
    return {'minLen':minLen, 'minValue': minValue, "allPaths": readyPath, 'alterMinLen':alterMinLen, 'alterMinValue':alterMinValue}

def rxInPath(path):
    return sum(path[1])

def pathLen(path):
    return len( path[0])

def pathValue(path):
    value=0
    for i in range(1, len(path[0])):
        value+= abs(path[0][i]-path[0][i-1])
    #print "PATH VALUE", value, "PATH", path
    return value


def growPath(pathOryg, allRxOryg):
    #print "PATH ORYG", pathOryg, "\t", 
    #path={path:[], nextCondition:3, rxToCheck=[ {}, {}, {} and all in further in allRxOryg]
    cond=pathOryg['nextCondition']
    pathWithThisIter=pathOryg['path'][:]
    if len(pathWithThisIter) > 200:
        print "\n ALL RX", allRxOryg, "\n"
        print pathWithThisIter, "\n"
        raise
    #
    if pathWithThisIter and cond == pathWithThisIter[-1]:
        print "==XXXX", pathWithThisIter, "CCC", cond, "ORYG", pathOryg
        raise
    pathWithThisIter.append(cond)
    #rxToCheckFromOryg=pathOryg['rxToCheck']
    rxToCheck=pathOryg['rxToCheck'][:]
    notMatchNearestRx=[]
    addedRx=[]
    itr=0
    changeAfterRx=[]
    allRxCopy=allRxOryg[:]
    if 'change' in pathOryg:
        changeAfterRx= pathOryg['change'][:]
        #print "XXXX", changeAfterRx, "ORYGPATH", pathOryg
        #raise
    while True:
        itr+=1
        newToCheck=[]
        test=[]
        for rx in rxToCheck:
            #print "RX", rxToCheck
            #if cond >= rx['condition'][0] and cond <= rx['condition'][1]:
            if cond in rx['condition']:
                if not rx in addedRx:
                    addedRx.append(rx)
                newToCheck.extend( [r for r in allRxCopy if r['product'] in rx['substrates'] ] )
            else:
                notMatchNearestRx.append(rx)
        #print "ADDED", len(newToCheck)
        if len(newToCheck) ==0:
            changeAfterRx.append( len(addedRx) )
            if len(notMatchNearestRx)> 0: # not full path
                #make new paths
                allowedValue=set()
                for newRx in notMatchNearestRx:
                    #print "COND", newRx['condition']
                    #for oneVal in range( newRx['condition'][0], 1+ newRx['condition'][1]):
                    for oneVal in newRx['condition']:
                        allowedValue.add(oneVal)
                pathListToRet=[]
                
                for c in allowedValue:
                    if c == cond: 
                        continue
                    #print "XXXX", pathWithThisIter, "CCC", c, "COND", cond
                    #modPaths={'path':path['path'].append( path['nextCondition']), 'nextCondition':None, 'rxToCheck':[]  ]
                    pathListToRet.append( ( {'path':pathWithThisIter, 'nextCondition':c, 'rxToCheck': notMatchNearestRx, 'change':changeAfterRx }, 'notReady' ) )
                #print "PATH LIST RET", pathListToRet
                return  pathListToRet
            else: 	#full path
                #print "FULL PATH", pathWithThisIter
                #print "XXXX", changeAfterRx
                #print "PPP", [ ( ( tuple(pathWithThisIter), "||", tuple(changeAfterRx) ), 'ready' )]
                return [ ( ( tuple(pathWithThisIter), tuple(changeAfterRx) ), 'ready' )]
            raise
        rxToCheck=newToCheck
        if itr > 500:
            print "\n path:", pathOryg, "ALL:", allRxOryg,
            #raise "ITR IS ABOVE 22500"
            return[ (False, 'ready')  ]
    ##ret [ (path,status), (path,status) ]



def getKgraph(dicData):

    allRx=set()
    for smi in dicData['results']:
        for rx in dicData['results'][smi]['history']:
            rxs= rx['smiles'].split('>>')
            subst= tuple(sorted(rxs[0].split('.')))
            prod=rxs[1]
            allRx.add( (subst, prod) )
        if 'alterHistories' in dicData['results'][smi] and dicData['results'][smi]['alterHistories']:
            for oneAlter in dicData['results'][smi]['alterHistories']:
                for rx in oneAlter:
                    rxs= rx['smiles'].split('>>')
                    subst= tuple(sorted(rxs[0].split('.')))
                    prod=rxs[1]
                    allRx.add( (subst, prod) )
    print "allRX", len(allRx)
    allSmi=dict()
    for (subst,prod) in allRx:
        for sub in subst:
            if not sub in allSmi:
                allSmi[sub]={'subst':0, 'prod':0}
            allSmi[sub]['subst']+=1
        if not prod in allSmi:
            allSmi[prod]={'subst':0, 'prod':0}
        allSmi[prod]['prod']+=1
    lifeCount=0
    noLifeSmi=[]
    for smi in allSmi:
        life="NO"
        if smi in lifeSmiles:
            life="LIFE"
            lifeCount+=1
            #print smi
        else:
            noLifeSmi.append( smi)
            import random
    numbers=set()
    while len(numbers) <27:
        rint=random.randint(0, len(noLifeSmi) )
        numbers.add(rint)
        print "LEN", len(numbers), len(noLifeSmi)
    for n in numbers:
        print noLifeSmi[n]
        #print life, '\t', smi, '\t', allSmi[smi]['subst'], '\t', allSmi[smi]['prod']
    return lifeCount




if __name__ == "__main__":
    import sys
    import cPickle as pickle
    dicData=pickle.load( open(sys.argv[1]) )
    makeCondVariant( dicData)
    #print getKgraph(dicData)

