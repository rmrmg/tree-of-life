### technical information
- all files are compresed with p7zip (https://www.7-zip.org/) 
- files are provided in two formats: json and python's pickle file (inside is one big dict)
- all pickle file was created under python2.7
- after uncompression files are 100 times bigger (check available space on your disk)
- loading pickle files into Python require a lot of RAM (e.g. 6 generations tree with all paths requires ~68GB of RAM!!!)

### data format
- data is store as big one python dict, which follwoing keys:
    - dict['status'] here is always just 'ok' string which confirm that the calculation successfully ended
    - dict['cycles'] list of cycles in graph (here cycles are truncated up to 5 reactions length)
    - dict['results'] dict where keys are smiles of reactants/products, e.g. dict['results']['C=C(C=NCC#N)c1ccccc1']. Such entry is also python dict with following keys:
        - parents - reactants which gave this product
        - otherParents - other reactans which can produce this compound
        - history - list of all reactions from initial reactants to this compounds (synthetic path). Each step in history is dict with following entries: 
            - smiles - reaction smiles
            - description - reaction name
            - reference - literature reference to the reaction type
            - additionalInfo - reaction's thermochemistry (dG, dS, dH calculated with Benson group contribution), and estimated (based on literature precedents) yield. Benson estimation is not very accurate for heterocycles so we switch to semi empirical calculation and dS, dG, dH may be 0 in the files.
            - condAB - reaction condition acidic, neutral, base or any combination of this
            - conditions - textual description of reaction condition
        - alterHistories - list of list of alternative synthetic path. Each step has the same format as in "history".
        - all other entries should be self explanatory
- data in json format correspond to python's pickle file 

### information about files
- tree_of_life_6generations_CNO tree with compounds up to 6th synthetic generation with mass limit 300g/mol and 5 reactans: water, amonia, nitrogen, methane and HCN. This tree has cycles and all alternative paths.
- tree_of_life_6generations_CNOS tree with compounds up to 6th synthetic generation with mass limit 300g/mol and 5 reactans: water, amonia, nitrogen, methane, HCN and H2S. This tree has cycles and all alternative paths. 
This part is splited into 3 parts with suffix _part_a _part_b and _part_c you can join it using cat command in unix, i.e.:  tree_of_life_6generations_CNOS.pickle.7z_part_a tree_of_life_6generations_CNOS.pickle.7z_part_b tree_of_life_6generations_CNOS.pickle.7z_part_c > tree_of_life_6generations_CNOS.pickle.7z



### hints
- if you are interesed in testing your hypothesis (e.g. what can be formed without some reactions etc, is it possible to make compound XXX without reaction YYY and/or reactant ZZZ and so on) the best way is to build graph from the data
- if you want to load data as graph, graph_tool seems to be the best choice (for example see gtCycleShortLenMultiproc.py)
- if you have any questions please contact with us