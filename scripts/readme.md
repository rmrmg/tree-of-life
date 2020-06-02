## example scripts which ilustrate how to works with tree-of-life data

## getSmilesAndMass.py 

the simplest example how to extract smiles and molecular weight

## nxCycleSearch.py 
load tree-of-life as graph using networkx library. The simples but not the most efficient way to work with graph.

## gtCycleShortLenMultiproc.py
 load tree-of-life as graph using graph_tool and search for cycles using multiprocesses. This is efficient was to search for cycles using all cores of your machine. 
Unfortunatelly installation of graph-tool can be more problematic (it cannot be done with pip) for details see https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions

## calcHeatOfFormation.py
calculate heat of formation with mopac. Mopac need to be installed in /opt/mopac (see http://openmopac.net/ for details how to get and install mopac)

## PlotProperties.py
Plots properties from CSV files deposited here. Assuming working directory is `scripts`, the following commands should reproduce figures from the paper.

* Figure 2a
```
 python  PlotProperties.py --descriptors heat_of_formation mass  Fig2a.png
```
* Figure 2b
```
python PlotProperties.py --descriptors LogP hb_acceptors Fig2b.png
```
* Figure S59
```
python PlotProperties.py FigS59.png
```
## CompareDistributions.py
Used for statistical tests (see `python CompareDistributions.py -h`). The following commands should reproduce results in TableS2 (note that Bootstrap results may differ by last digit due to the stochastic character of the method).

* Number of condition changes
```
python CompareDistributions.py --descriptors conditions
```

* LogP
```
python CompareDistributions.py --descriptors LogP 
```

* 2-D distributions
   * LogP vs Nuber of Hydrogen Bond Donors
```
python CompareDistributions.py --descriptors LogP hb_donors --two_d
```

   * LogP vs Number of Rings
```
python CompareDistributions.py --descriptors LogP rings --two_d
```
   * Hydrogen Bond Donors vs Hydrogen Bond Acceptors
```
python CompareDistributions.py --descriptors hb_donors hb_acceptors --two_d
```
