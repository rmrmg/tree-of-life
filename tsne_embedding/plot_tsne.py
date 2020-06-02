#========= Arguments
import argparse
parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--biotic', type=str, default='../feature_elimination/biotic_rdkit.npz', help=' ')
parser.add_argument('--abiotic', type=str, default='../feature_elimination/abiotic_rdkit.npz', help=' ')
parser.add_argument('--feature_names', type=str, default='../feature_elimination/rdkit_desc_names_sorted.txt', help=' ')
parser.add_argument('--selection', type=str, default='selected_features.txt', help='List of selected features')
parser.add_argument('output', type=str, help='Name of the output file')
args=parser.parse_args()
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s:%(levelname)s: %(message)s')
#========= AXES DATA
logging.info('Start')

with open(args.feature_names, 'r') as f:
    desc_names = [x.strip() for x in f]
 
if not (args.selection is None): 
    with open(args.selection, 'r') as f:
        selection = [x.strip() for x in f]
    desc_idx = [i for i,x in enumerate(desc_names) if x in selection]
else:
    desc_idx = [i for i in range(len(desc_names))]
 
#========= OTHER
import numpy as np
import gzip
import pickle as pic
#from matplotlib import use
#use('Agg')
from matplotlib import pyplot as plt
from sklearn import manifold
seed = np.random.RandomState(seed=3)

#======== READ ARRAYS
with open('clusters.pic','rb') as f: clusters = pic.load(f)

with gzip.open(args.biotic,'rb') as f:biotic_rdkit=np.load(f)[:,desc_idx]
with gzip.open(args.abiotic,'rb') as f:abiotic_rdkit=np.load(f)

abiotic_rdkit = np.array([abiotic_rdkit[x[0],desc_idx] for x in clusters])

sizes=(np.loadtxt('cluster_sizes.txt')/2).astype(int)

Nabiotic = abiotic_rdkit.shape[0]
X=np.vstack([abiotic_rdkit, biotic_rdkit])

std=X.std(axis=0)
u=X.mean(axis=0)
std=np.where(std>0,std,1)
X=(X-u)/std

#tsne
logging.info('Data read')
tsne=manifold.TSNE(n_components=2, metric = 'euclidean')
newX = tsne.fit_transform(X)
logging.info('TSNE Done')
np.save('%s-tsne.npy'%args.output[:-4], newX)
logging.info('TSNE written')
#======= Plot
fig, ax = plt.subplots()
ax.scatter(newX[:Nabiotic, 0], newX[:Nabiotic, 1], color='#D6DCE5', s=sizes, lw=0, label='Abiotic molecules',zorder=0)
handle_abiotic = ax.scatter([],[],color='#D6DCE5', s=100)
ax.scatter(newX[Nabiotic:, 0], newX[Nabiotic:, 1], color='#2F5597',marker='*', s=50, lw=0, label='Biotic molecules',zorder=2)
handle_biotic = ax.scatter([],[], color='#2F5597',marker='*',s=100)
#ax.set_xlabel(labels[0],fontsize=12)
#ax.set_ylabel(labels[1],fontsize=12)
ax.legend(scatterpoints=1, loc='best', shadow=False, handles=[handle_biotic, handle_abiotic], labels=['Biotic','Abiotic'])
plt.savefig(args.output,dpi=300)
logging.info('Plot Done')
