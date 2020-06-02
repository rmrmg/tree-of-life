from scipy.special import kolmogorov
from scipy.stats import chi2 as scipy_chi2
from math import sqrt
from random import randint
from numpy.random import multivariate_normal
from numpy import array,zeros,ones,where
import numpy as np
from rdkit.Chem import CanonSmiles, MolFromSmiles
from rdkit.Chem.Descriptors import RingCount
from multiprocessing import Pool

def canonize(sml):
   return CanonSmiles(CanonSmiles(sml))

def chi2(x,v): return 1-scipy_chi2.cdf(x,v)

def chi2_stat(X,Y):
    N_X=float(len(X))
    N_Y=float(len(Y))
    kx=sqrt(N_Y/N_X)
    ky=sqrt(N_X/N_Y)
    R=max( max(X),max(Y) )
    S=0
    for i in range(R+1):
        #get counts back
        nx=len([x for x in X if x==i])
        ny=len([x for x in Y if x==i])
        if nx+ny>0:
            S+=(kx*nx-ky*ny)**2/(nx+ny)
    return S,R

def make_cum(X):
    R=max(X)
    N=float(len(X))
    return [len([x for x in X if x<=i])/N for i in range(R+1)],N

def ks_test(X,Y):
    cX,Nx=make_cum(X)
    cY,Ny=make_cum(Y)
    S=0
    lX=len(cX)
    lY=len(cY)
    
    for i in range(max(lX,lY)):
        if i==0:
            d=abs(cX[i]-cY[i])
        elif i<lX and i<lY:
            d=max(abs(cX[i]-cY[i]), abs(cX[i-1]-cY[i]), abs(cX[i]-cY[i-1]))
        elif i<lX:d=max(abs(cX[i-1]-1.0),abs(cX[i]-1.0))
        elif i<lY:d=max(abs(cY[i-1]-1.0),abs(cY[i]-1.0))
        else: d=0
        
        if d>S:S=d
    return S*sqrt(Nx*Ny/(Nx+Ny))


def bootstrap_internal(args):
    NX, NY, ref, stat, combined = args
    newX=[combined[randint(0,NX+NY-1)] for x in range(NX)]
    newY=[combined[randint(0,NX+NY-1)] for x in range(NY)]
    V=stat(newX,newY)
    if type(V).__name__=='tuple':V=V[0]
    if V>=ref: 
      return 1.0
    return 0.0


def bootstrap(X,Y,nb=1000,stat=lambda x:0.001, threads=1):
    combined = np.hstack([X,Y])
    NX=len(X) 
    NY=len(Y)
    pval=0
    ref=stat(Y,X)
    if type(ref).__name__=='tuple':ref=ref[0]

    if threads==1:
      internal = lambda x: bootstrap_internal((NX, NY, ref, stat, combined))
      pval = sum([internal(_) for _ in range(nb)])
    else:
      p=Pool(threads)
      data_gen = ((NX, NY, ref, stat, combined) for _ in range(nb))
      pval = sum(p.map(bootstrap_internal, data_gen))
    return pval/nb
    
def mean(X):return sum(X)/float(len(X))

def _make_covariance(X):
    R=len(X)
    x1=X[:-1]
    #x2=x[1:]
    z1=zeros((R-1,R-1))
    #z2=zeros((R-1,R-1))
    for i in range(R-1):
        for j in range(i,R-1):
            z1[i,j]=min(x1[i],x1[j])-x1[i]*x1[j]
            #z2[i,j]=min(x2[i],x2[j])-x2[i]*x2[j]
            z1[j,i]=z1[i,j]
            #z2[j,i]=z2[i,j]
    return z1#,z2

def multi_KS(nX,nY,n=100):
    X,x=make_cum(nX)
    Y,y=make_cum(nY)
    #X=X[2:]
    #Y=Y[2:]
    print( X,Y)
    Z1=_make_covariance(X)
    Z2=_make_covariance(Y)
    #print('COV DONE')
    print(Z1)
    #print(Z2)
    n1=len(X)-1
    n2=len(Y)-1
    U1=zeros(n1)
    U2=zeros(n2)
    KS=ks_test(nX,nY)
    o1=ones(len(X)-1)*KS
    print( o1)
    o2=ones(len(Y)-1)*KS
    z1=multivariate_normal(U1,Z1,size=n,check_valid='warn')-o1
    print(max(z1[0]))
    z2=multivariate_normal(U2,Z2,size=n,check_valid='warn')-o2
    z1=where(z1<0,1,0).sum(axis=1)
    print(z1.shape)
    z2=where(z2<0,1,0).sum(axis=1)
    z1=where(z1==n1,1,0).mean()
    z2=where(z2==n2,1,0).mean()
    return z1,z2
    

def counts_to_bin_ids_representing_data(counts):
   result=[]
   for i, count in enumerate(counts):
      result.extend([i for _ in range(count)])
   return result


if __name__=='__main__':
   import argparse
   import pandas as pd
   parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('--properties', type=str, 
         default='../molecule_data/selected_properties.csv',
         help='path to selected_properties.csv')
   parser.add_argument('--conditions', type=str,
         default='../molecule_data/minimalCondtionsChanges.csv', 
         help='path to minimalCondtionsChanges.csv')
   parser.add_argument('--nbootstrap',type=int,default=1000, help='Number of bootstrap samples')
   parser.add_argument('--nthreads',type=int,default=1, help='Number of threads')
   parser.add_argument('--two_d',action='store_true', help='Compare 2D distributions (2 descriptors must be provided)')
   parser.add_argument('--descriptors', type=str, nargs='+', 
         default=['LogP'],
         choices=['hb_donors', 'hb_acceptors', 'LogP', 'mass', 'conditions', 'rings'],
         help='Properties to compare')
   args=parser.parse_args()
   
   keys= {'mass':'ExactMolWt','LogP':'MolLogP','hb_acceptors':'NumHAcceptors',
          'hb_donors':'NumHDonors', 'heat_of_formation':'heat', 
           'conditions':'conditions', 'rings':'rings'}
   properties = pd.read_csv(args.properties, sep=';')
   #properties.smiles = properties.smiles.apply(canonize)
   properties.sort_values('smiles', inplace=True, ignore_index=True )

   if 'rings' in args.descriptors:
      properties = properties.assign(rings=properties.smiles.apply(lambda x: RingCount(MolFromSmiles(x))))

   if 'conditions' in args.descriptors:
      conditions = pd.read_csv(args.conditions, sep=';')
      #conditions.smiles = conditions.smiles.apply(canonize)
      conditions.sort_values('smiles',inplace=True, ignore_index=True )
      conditions = conditions[conditions.smiles.isin(properties.smiles)]

      properties = properties[properties.smiles.isin(conditions.smiles)]
      properties = properties.assign(conditions=conditions[conditions.smiles.isin(properties.smiles)]['minimal number of condition changes '])

   nb = args.nbootstrap

   biotic_series, abiotic_series = [], []
   do_2d = (args.two_d and len(args.descriptors)==2)
   for desc in args.descriptors:
      biotic_data = properties[properties.biotic_flag==1][keys[desc]]
      biotic_series.append(biotic_data)

      abiotic_data = properties[properties.biotic_flag==0][keys[desc]]
      abiotic_series.append(abiotic_data)
      
      if not do_2d:
         print('Property :%s'%desc)
         bins = 10 if desc in ['mass','LogP'] else np.arange(properties[keys[desc]].max()+1)
         _, bins = np.histogram(properties[keys[desc]], bins=bins)
      
         biotic_data, _ = np.histogram(biotic_data, bins=bins)
         
         #now let's make a population of bin ids - each id should repeat the same number of time 
         #as count within bin
         #thus, during sampling, given bin will be represented proprotionally to its occupation
         #it is effectively the same as assining each record with the corresponding bin ID
         biotic_data = counts_to_bin_ids_representing_data(biotic_data)
   
         abiotic_data, _ = np.histogram(abiotic_data, bins=bins)
         abiotic_data = counts_to_bin_ids_representing_data(abiotic_data)
   
         #Chi2 p-values      
         c,v=chi2_stat(biotic_data, abiotic_data)
         cp=chi2(c,v)
         print('   Chi2: ',cp)
         
         #KS p-values
         k=ks_test(biotic_data, abiotic_data)
         kp=kolmogorov(k)
         print('   KS: ',kp)
         print('\n   Bootstrap')
          
         print('   KS: ',bootstrap(biotic_data, abiotic_data, nb=nb, stat=ks_test, threads=args.nthreads))
         print('   Chi2: ',bootstrap(biotic_data, abiotic_data, nb=nb, stat=chi2_stat, threads=args.nthreads))
   if do_2d:
      print('2D: %s-%s'%tuple(args.descriptors)) 
      _, binx, biny = np.histogram2d(properties[keys[args.descriptors[0]]], properties[keys[args.descriptors[1]]])
      bins=(binx, biny)
      biotic_data, _, _ = np.histogram2d(*biotic_series, bins=bins)
      biotic_data = counts_to_bin_ids_representing_data(biotic_data.reshape(-1).astype(int))

      abiotic_data, _, _ = np.histogram2d(*abiotic_series, bins=bins)
      abiotic_data = counts_to_bin_ids_representing_data(abiotic_data.reshape(-1).astype(int))
   
      #Chi2 p-values      
      c,v=chi2_stat(biotic_data, abiotic_data)
      cp=chi2(c,v)
      print('   Chi2: ',cp)
         
      print('\n   Bootstrap')
      print('   Chi2: ',bootstrap(biotic_data, abiotic_data, nb=nb, stat=chi2_stat, threads=args.nthreads))

      
      
