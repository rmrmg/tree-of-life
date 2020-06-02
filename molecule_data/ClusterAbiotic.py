from rdkit import DataStructs, Chem
import numpy as np
from rdkit.ML.Cluster import Butina
from rdkit.Chem import AllChem
import pickle as pic
import pandas as pd
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s:%(message)s')


data = pd.read_csv('selected_properties.csv', sep=';')
non_life_mols= data[data.biotic_flag==0].smiles.values

logging.info('data read')
fgps=[AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(mol),2,useFeatures=1,nBits=4096) for mol in non_life_mols]


def ClusterFps(fps,cutoff=0.5):
    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    print_mem(stage='Dist & Clust usage')
    del dists
    return cs

logging.info('fgp done')

cs=ClusterFps(fgps)
with open('clusters_test.pic','wb') as f: pic.dump(cs,f)

logging.info('clusters written')

is_centroid = np.zeros(len(non_life_mols)).astype(int)
clust_idx = np.zeros(len(non_life_mols)).astype(int)


for idx,C in enumerate(cs):
   is_centroid[C[0]]=1
   clust_idx[list(C)]=idx


df = pd.DataFrame(dict(clust_id=clust_idx, is_centroid=is_centroid, smiles=non_life_mols))
df=df.assign(temp=df.clust_id*10-df.is_centroid)
df.sort_values('temp', inplace=True)
df = df.drop('temp', 1)

df.to_csv('abiotic_clusters.csv', sep=';', index=False)
