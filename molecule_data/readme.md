# Information about files
All files are in CSV format (with semicolon as a  field separator).

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

## t-SNE Embeding of molecules in Fig S60
Note that due to a stochastic character of the algorithm, actual coordiantes may change when computed the second time.
* file: `fig_s60_crd.csv`
* columns: `biotic_flag`, `cluster_id` and `smiles` have the same meaning as above
* only cluster centroids are presented
* `cluster_size` denotes the size of the corresponding cluster
* `x` and `y` denote coordinates in the embedded space
* for biotic compounds, `cluster_id` and `cluster_size` does not apply, thus are assigned with values `not apply`

## Path statistics

`KinKout.csv`  contains ***in degree*** (how many reactions leads to particular molecule) and ***out degree*** (number of reactions in which molecules acts as substrate ) of nodes (molecules) in graph (tree of life).

