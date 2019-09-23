import sys
import cPickle as cp

#the simplest example how to work with trees of life
#here you iterate over all compounds and print its smiles and molecular mass

dane=cp.load( open(sys.argv[1]) )
for i in dane['results']:
    print i, "has mass", dane['results'][i]['mass']
