#!/usr/bin/python2
from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem
import sys, os
from multiprocessing import Pool

th1=0.8
th2=0.4

class loopEscape(Exception):
    pass

def retxyz(smiles, fftype="MMFF"):
	mol=Chem.MolFromSmiles(smiles)
	mol=Chem.AddHs(mol)
	confs=500
	numThreads=12
	RMSThreshold=th1
	useBK=False #True #use basic knowledge about
	useExpTors=False #True #use experimental value for tortian angles
	howmany=Chem.rdDistGeom.EmbedMultipleConfs(mol, numConfs=confs, randomSeed=-1,  numThreads=numThreads, 
	  pruneRmsThresh=RMSThreshold, useBasicKnowledge=useBK, useExpTorsionAnglePrefs=useExpTors )
	#print(len(howmany))
	minimE=None
	minimConf=None
	allE=dict()
	for confid in range(len(howmany)):
		if fftype=="UFF":
			result= AllChem.UFFOptimizeMolecule(mol, maxIters=9500, vdwThresh=80.0, confId=confid, ignoreInterfragInteractions=False) 
			ff = AllChem.UFFGetMoleculeForceField(mol, confId=confid)
		elif fftype=="MMFF":
			result= AllChem.MMFFOptimizeMolecule(mol, maxIters=12540, confId=confid, ignoreInterfragInteractions=False)
			mmffProp=AllChem.MMFFGetMoleculeProperties(mol, mmffVerbosity=0)
			ff = AllChem.MMFFGetMoleculeForceField(mol, mmffProp, confId=confid, nonBondedThresh=80)
			#energy_value = ff_min.CalcEnergy()
		energy = ff.CalcEnergy()
		#print("energy of this conf: "+str(energy) )
		if minimE == None or minimE > energy:
			minimE=energy
			minimConf=confid
		allE[confid]=energy
	return (mol, allE)

def uniqConformers(mol, allE):
    edict=dict()
    #minimalE = min(allE.values() )
    for confID in allE:
        ener=round(allE[confID],3)
        if not ener in edict:
            edict[ener]=[]
        edict[ener].append(confID)
    #print edict
    symbols=[ a.GetSymbol() for a in mol.GetAtoms() ]
    tmp=tuple(set(symbols))
    symbolIdxs=dict()
    for i,symb in enumerate(symbols):
        if not(symb in symbolIdxs):
            symbolIdxs[symb]=[]
        symbolIdxs[symb].append(i)
    atomGroups=[]
    for i in symbolIdxs:
        atomGroups.append( symbolIdxs[i] )
    #print symbols, symbolIdxs
    #print atomGroups
    #identic=dict()
    identic=[]
    for oneE in edict:
        conformerIDs=edict[oneE]
        for idx1 in range( len(conformerIDs)):
            confId1=conformerIDs[idx1]
            conf1=mol.GetConformer(confId1)
            for idx2 in range(idx1+1, len(conformerIDs)):
                confId2=conformerIDs[idx2]
                conf2=mol.GetConformer(confId2)
                try:
                    for gr1id in range(len(atomGroups)):
                        for gr2id in range(len(atomGroups)):
                            if areTheSame(conf1, atomGroups[gr1id], atomGroups[gr2id], conf2, atomGroups[gr1id], atomGroups[gr2id] ):
                                #identic.append( (confId1, confId2) )
                                #if not (confId1 in identic):
                                #    identic[confId1]=[]
                                #identic[confId1].append(confId2)
                                identic.append(confId2)
                                raise loopEscape
                except loopEscape:
                    #print "identical", confId1, confId2
                    pass
    uniq=[]
    for confId in allE:
        if not confId in identic:
            uniq.append(confId)
    #for
    return uniq

def areTheSame(conf1, c1ag1, c1ag2, conf2, c2ag1, c2ag2):
    conf1dist=[]
    conf2dist=[]
    epsilon=0.01

    for atm1id in c1ag1:
        for atm2id in c1ag2:
            if atm1id == atm2id:
                continue
            dist=conf1.GetAtomPosition(atm1id).Distance(conf1.GetAtomPosition(atm2id) )
            conf1dist.append( dist)
    conf1dist.sort()
    for atm1id in c2ag1:
        for atm2id in c2ag2:
            if atm1id == atm2id:
                continue
            dist=conf2.GetAtomPosition(atm1id).Distance(conf2.GetAtomPosition(atm2id) )
            conf2dist.append( dist)
    conf2dist.sort()
    ##print  len(conf1dist)
    for idx in range( len(conf1dist) ):
        ##print conf1dist[idx], conf2dist[idx]
        if abs(conf1dist[idx]-conf2dist[idx]) > epsilon:
            return False
    return True

def saveXyz(mol, smiles, fftype, confid, E, doSave=True):
	if doSave:
		FILE=open('res'+smiles+'_'+fftype+'.xyz', 'a')
	molstr=Chem.MolToMolBlock(mol, confId=confid ).split('\n')
	#print 'molstr', molstr
	#print '---', molstr[3]
	ifpr=False
	atmnum=999 #anything more that 3 where true value is located
	#allConformers=[]
	newConf=[]
	for lineid in range(len(molstr)):
		if lineid == 3:
			#print 'lineid', molstr[lineid] 
			atmnum= int( molstr[lineid].split()[0])
			if doSave:
				print >>FILE, atmnum
				print >>FILE, E
			ifpr=True
		elif ifpr:
			line=molstr[lineid].split()
			#print 'LINE::', line
			topr=line[3]+ " " + line[0] + " " + line[1] +" "+ line[2]
			newConf.append( topr.strip() )
			if doSave:
				print >>FILE, topr
		if lineid >= 3+atmnum:
			#allConformers.append(newConf)
			break
	return newConf



#header="AM1 OPT FREQ\n20.2658\n"
def makeMopacCalc(geomConfId, theory='PM6-D3H4X', prefix='inp', doOpt=True):
    #print "XX", geomConfId
    confNum, geom = geomConfId
    #print("CC", confNum)
    if doOpt:
        headerLines=theory #+ " NOOPT " #+' OPT '
    else:
        headerLines=theory+ " NOOPT " 

    fw = open('tmp_%i.mop' % confNum, 'w')
    #print >>fw, headerLines
    #print >>fw, "comment" # thisConformer[1][:-1]
    #print >>fw, ''
    print(headerLines, file=fw)
    print( "comment", file=fw)
    print('', file=fw)
    for ln in geom:
        #print >>fw, ln
        print(ln, file=fw)
    fw.close()
    #print "START", confNum
    os.system('LD_LIBRARY_PATH=/opt/mopac /opt/mopac/MOPAC2016.exe tmp_%i.mop 2> /dev/null' % confNum)
    #print "END"
    energy = [l for l in open('tmp_%i.out' % confNum) if 'TOTAL ENERGY' in l][0].split()[-2] #eV
    formationEnergy = [l for l in open('tmp_%i.out' % confNum) if 'FINAL HEAT OF FORMATION ' in l][0].split()[5] #kcal/mol
    #print "E", float(energy)
    return {'totalE_eV': float(energy), 'heatForm_kcal/mol':float(formationEnergy) }


def getMopacEnergy(smiles, fftype ):
	res=retxyz(smiles, fftype=fftype)
	confId=uniqConformers(res[0], res[1])
	#confId=res[1].keys()
	#print res[1]
	minE=min( res[1].values() )
	allConf=[]
	for idx in confId:
		E=round(res[1][idx] - minE,4)
		if E > 5:
			continue
		allConf.append( saveXyz(res[0], smiles, fftype, idx, E, doSave=False) )
	#p=Pool(2)
	minimalE = min([ hf['heatForm_kcal/mol'] for hf in map(makeMopacCalc, enumerate(allConf) )])
	os.system('rm tmp*')
	return minimalE

if __name__ == "__main__":
	smilesFile = open(sys.argv[1])
	fftype='MMFF'
	#if len(sys.argv) ==3:
	#	fftype=sys.argv[2]
	for smiles in smilesFile:
		try:
			print(smiles.strip(), getMopacEnergy( smiles, fftype))
		except:
			print("CANNOT CALC", smiles.strip())


