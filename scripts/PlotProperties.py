#========= Arguments
import argparse
parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--properties', type=str, 
      default='../molecule_data/selected_properties.csv',
      help='path to selected_properties.csv')
parser.add_argument('--heat_of_formation', type=str,
      default='../molecule_data/heatOfFormation.csv', 
      help='path to heatOfFormation.csv')
parser.add_argument('--abiotic_clusters', type=str, 
      default='../molecule_data/abiotic_clusters.csv',
      help='path to abiotic_clusters.csv')
parser.add_argument('--use_all_abiotic', action='store_true',
      help='Use all abiotic compounds or just cluster centroids?')
parser.add_argument('--descriptors', type=str, nargs=2, 
      default=['hb_donors', 'hb_acceptors'],
      choices=['heat_of_formation', 'hb_donors', 'hb_acceptors', 'LogP', 'mass'],
      help='Descriptors to plot in the order X,Y (use names of columns in CSV files')
parser.add_argument('--alpha', type=float, default=1., 
      help='Alpha value for abiotic compounds (1-opaque, 0-totally transparent)')
parser.add_argument('output', type=str, help='path to the output file')
args=parser.parse_args()

import pandas as pd
properties = pd.read_csv(args.properties, sep=';').sort_values('smiles')
heat = pd.read_csv(args.heat_of_formation, sep=';').sort_values('smiles')
heat = heat[heat.smiles.isin(properties.smiles)]
heat.drop_duplicates('smiles',inplace=True)
properties = properties.assign(heat=heat[' heat of formation [kcal/mol]'].values)
properties = properties.assign(size=[10 for _ in properties.smiles])
abiotic_clusters = pd.read_csv(args.abiotic_clusters, sep=';').sort_values('smiles')

labels = {'hb_donors':'# Hydrogen Bond Donors', 'hb_acceptors':'# Hydrogen Bond Acceptors', 
          'LogP':'logP', 'mass': 'Molecular mass', 'heat_of_formation':'Heat of formation [kcal/mol]'}
keys= {'mass':'ExactMolWt','LogP':'MolLogP','hb_acceptors':'NumHAcceptors', 'hb_donors':'NumHDonors', 'heat_of_formation':'heat'}

Xk, Yk = [keys[x] for x in args.descriptors]
Xl, Yl = [labels[x] for x in args.descriptors]

abiotic_c = 'k'

#limit data to cluster_centroids
if not args.use_all_abiotic:
   part = abiotic_clusters[abiotic_clusters.is_centroid==1]
   sizes = [(abiotic_clusters.cluster_id==x).sum()/2 for x in part.cluster_id]
   properties = properties[(properties.smiles.isin(part.smiles)) | (properties.biotic_flag==1)]
   properties.loc[properties.biotic_flag==0, 'size'] = sizes
   abiotic_c = '#D6DCE5'

 
#========= OTHER
import numpy as np
import gzip
import pickle as pic
#from matplotlib import use
#use('Agg')
from matplotlib import pyplot as plt
seed = np.random.RandomState(seed=3)

#======= Plot
fig, ax = plt.subplots()
ax.scatter(properties[properties.biotic_flag==0][Xk],
           properties[properties.biotic_flag==0][Yk], 
           color=abiotic_c, 
           alpha=args.alpha,
           s=properties[properties.biotic_flag==0]['size'], 
           lw=0, label='Abiotic molecules', zorder=0)
handle_abiotic = ax.scatter([],[],color='#D6DCE5', s=100)
ax.scatter(properties[properties.biotic_flag==1][Xk],
           properties[properties.biotic_flag==1][Yk], 
           color='#2F5597',
           marker='*', 
           s=50, 
           lw=0, label='Biotic molecules', zorder=2)
handle_biotic = ax.scatter([],[], color='#2F5597',marker='*',s=100)
ax.set_xlabel(Xl,fontsize=12)
ax.set_ylabel(Yl,fontsize=12)
ax.legend(scatterpoints=1, loc='best', shadow=False, handles=[handle_biotic, handle_abiotic], labels=['Biotic','Abiotic'])
plt.savefig(args.output,dpi=300)
