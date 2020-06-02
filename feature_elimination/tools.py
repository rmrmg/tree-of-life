import numpy as np
import gzip
import pickle as pic

def pickle_gz(name, obj):
    with gzip.open(name, 'wb') as f:
        pic.dump(obj, f)


def unpickle_gz(name):
    with gzip.open(name, 'rb') as f:
        return pic.load(f)


def load_npz(name):
    with gzip.open(name, 'rb') as f:
        data = np.nan_to_num(np.load(f))
        data=np.where(data<1e9,data,1e9)
        return data
    
    
def load_txt(name):
    with open(name, 'r') as f:
        return np.array(f.readlines())


def save_npz(name,data):
    with gzip.open(name, 'wb') as f:
        np.save(f,data)


def save_txt(name, data):
    np.savetxt(name, data, fmt='%s')
    
    
def load(name):
    name_list = name.split('.')
    if name_list[-1]=='npz':
        return load_npz(name)
    return load_txt(name)

def _flatten(array):
    N=1
    for x in array.shape:
        N*=x
    return array.reshape(N)


def clean_zeros(array1D, forbidden=[],epsilon=1e-6, renormalize=False):
   for i in range(len(array1D)):
      if i in forbidden: continue
      if array1D[i]==0:
         array1D[i]=epsilon
   if renormalize:
      array1D/=array1D.sum()
   return array1D


def find_names_correlated_with_Wt(X, names, th, key='Wt'):
    ref_indices = [i for i,x in enumerate(names) if key in x]
    result = []
    if (not (th is None)) and (ref_indices!=[]):
      idx = ref_indices[0]
      for I, name in enumerate(names):
         R = abs(np.corrcoef(X[:,idx], X[:,I])[0,1])
         if R>th:
            result.append(name)
    
    return result


def rdkit_eliminate_features(X, names, forbidden_names=['Wt'], corr_to_mass_th=None):
    #1 eliminate zero std
    std = X.std(axis=0)
    idx = np.where(std>0)[0]
    #2 eliminate weight and other stuff
    forbidden_names+= find_names_correlated_with_Wt(X, names, corr_to_mass_th)   
    if forbidden_names!=[]:
        idx = [i for i in idx if all([x not in names[i] for x in forbidden_names])]
    return X[:,idx], [names[i] for i in idx]


def calc_rates(clf, X, Y):
    pred=clf.predict(X)
    pos_n = Y.sum()
    neg_n = len(Y)-pos_n
    true_pos = pred.dot(Y)
    true_neg = (1-pred).dot(1-Y)
    return true_pos/pos_n, true_neg/neg_n
