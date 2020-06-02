# Information about files
All files are in CSV format (with semicolon as a  field separator).

Packages required to launch scripts:
* Scipy (version >=1.4.1)
* Rdkit (version >=2019.03.4)
* Pandas (version >=1.0.2)
* Matplotlib (version >=3.1.3)

## Statistics of generated compouds

* `selected_properties.csv` - computed properties for both biotic and abiotic compounds
 columns:
    - `smiles`: SMILES string of molecule
    - `biotic_flag`: `1`-biotic compound, `0`-abiotic compound
    -  `ExactMolWt`: Molecular mass (atomic units)
    -  `MolLogP`: LogP calculated with Crippen method
    -  `NumHAcceptors`: number of hydrogen bond acceptors
    -  `NumHDonors`: number of hydrogen bond donors
* `abiotic_clusters.csv` - clusterisation of abiotic compounds 
 columns:
 -`cluster_id`: ID of cluster
 -`is_centroid`: `1` for cluster centroids, `0` for other members of the cluster
 -`smiles`: SMILES of the molecule
Clusters are sorted according to their size (from the largest to smallest), with cluster centroid listed first.
The clusterisation result can be reproduced with `ClusterAbiotic.py` script (`selected_properties.csv` should be present in working directory). Note that execution of this scipt require about 21GB of memory.

## t-SNE Embeding of molecules in Fig S60
Note that due to a stochastic character of the algorithm, actual coordiantes may change when computed the second time.
* file: `fig_s60_crd.csv`
* columns: `biotic_flag`, `cluster_id` and `smiles` have the same meaning as above
* only cluster centroids are presented
* `cluster_size` denotes the size of the corresponding cluster
* `x` and `y` denote coordinates in the embedded space
* for biotic compounds, `cluster_id` and `cluster_size` does not apply, thus are assigned with values `not apply`


## heat of formations

`heatOfFormation.csv` contains heat of formation (in kcal/mol) calculated using PM6-D3H4X semiempirical method

## condition changes along pathway

`minimalCondtionsChanges.csv` contains information about minimal changes reaction condition, id est in which condition (from strongly acidic to strongly basic) reactions 
should occur to minimalize number of condition change
