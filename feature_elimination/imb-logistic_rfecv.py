import argparse
from tools import load_npz, np, save_npz, rdkit_eliminate_features, calc_rates
parser=argparse.ArgumentParser()
parser.add_argument('--biotic', type=str, default='biotic_rdkit.npz')
parser.add_argument('--abiotic', type=str, default='abiotic_rdkit.npz')
parser.add_argument('--corr_th', type=float, default=None)
parser.add_argument('--standarize', action='store_true')
parser.add_argument('--exclude', type=str, nargs='+', default=[])
parser.add_argument('--feature_names', type=str, default='rdkit_desc_names_sorted.txt')
parser.add_argument('--log', type=str, default=None)
args=parser.parse_args()
from sklearn.metrics import balanced_accuracy_score, make_scorer
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.model_selection import cross_val_score
from sys import stderr
import logging
#logging.basicConfig(level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')

lvl=logging.INFO
formatter = logging.Formatter('%(asctime)s : %(levelname)s : %(message)s')
logger = logging.RootLogger(logging.INFO)

stderr_handler = logging.StreamHandler(stderr)
stderr_handler.setLevel(logging.INFO)
stderr_handler.setFormatter(formatter)
stderr_handler.setLevel(lvl)
logger.addHandler(stderr_handler)


# initial configuration
SYM_NORM = True
if not(args.log is None):
   fh = logging.FileHandler(filename=args.log)
   fh.setLevel(logging.INFO)
   fh.setFormatter(formatter)
   logger.addHandler(fh)



def calc_rates(clf, X, Y):
    pred=clf.predict(X)
    pos_n = Y.sum()
    neg_n = len(Y)-pos_n
    true_pos = pred.dot(Y)
    true_neg = (1-pred).dot(1-Y)
    return true_pos/pos_n, true_neg/neg_n


with open(args.feature_names, 'r') as f:
    names = [x.strip() for x in f.readlines()]

biotic = load_npz(args.biotic)
abiotic = load_npz(args.abiotic)
Y = np.hstack([np.ones(biotic.shape[0]), np.zeros(abiotic.shape[0])])
X = np.vstack([biotic, abiotic])

X, names = rdkit_eliminate_features(X, names, forbidden_names = ['Wt']+args.exclude, corr_to_mass_th=args.corr_th)

standarizer = StandardScaler()
standarizer.fit(X)
if args.standarize:
    X=standarizer.transform(X)

scorer=make_scorer(balanced_accuracy_score)
log = LogisticRegression(penalty='l2', class_weight='balanced', solver='lbfgs', C=1.0, max_iter=1000)

Nfeatures = X.shape[1]
current_idx = np.arange(Nfeatures)
best_features = {'idx':[], 'av_score':0.0, 'std_score':0.0}
best_and_smallest = {'idx':[], 'av_score':0.0, 'std_score':0.0}

train_idx, test_idx = train_test_split(np.arange(X.shape[0]), test_size=0.2, stratify=Y)

logger.info('Commencing elimination')
for _ in range(Nfeatures, 1, -1):
        score = cross_val_score(log, X[:,current_idx], Y, cv=5, scoring=scorer)
        log2=LogisticRegression(penalty='l2', class_weight='balanced', solver='lbfgs', C=1.0, max_iter=1000)
        log2.fit(X[train_idx,:][:,current_idx],Y[train_idx])
        TP, TN = calc_rates(log2, X[test_idx,:][:,current_idx], Y[test_idx])
        av, std = score.mean(), score.std()/2
        if av > best_features['av_score']:
            best_features['idx']=current_idx[:]
            best_features['av_score']=av
            best_features['std_score']=std
        elif av>= (best_features['av_score']-best_features['std_score']):
            best_and_smallest['idx']=current_idx[:]
            best_and_smallest['av_score']=av
            best_and_smallest['std_score']=std
        logger.info('Nfeat :  %5i   bacc: %8.3f+- %8.3f  TP %8.3f  TN %8.3f  (best at %i : %8.3f+- %8.3f)'%(\
            len(current_idx), av, std, TP, TN, len(best_features['idx']), best_features['av_score'], best_features['std_score']))
        #update idx:
        log.fit(X[:,current_idx],Y)
        coef = abs(log.coef_.reshape(-1))
        if len(coef)==2: 
            break
        min_coeff = min(coef)
        current_idx = [x for i,x in enumerate(current_idx) if coef[i] > min_coeff]

logger.info('Best: %i, %8.3f+- %8.3f'%(len(best_features['idx']), best_features['av_score'], best_features['std_score']))
                                                                                                            
logger.info('Best and smallest: %i, %8.3f+- %8.3f'%(len(best_and_smallest['idx']), best_and_smallest['av_score'], best_and_smallest['std_score']))

bas_idx = best_and_smallest['idx']                                                                                                                                     
log.fit(X[:,bas_idx], Y)
names_f = list(zip([names[i] for i in bas_idx], log.coef_.reshape(-1)))
names_f.sort(key=lambda x:x[1], reverse=1)
logger.info('intercept: '+str(log.intercept_))
logger.info('selected features:')
for x,v in names_f:
    logger.info('%20s %10s'%(str(x),str(v)))
                                                                                                            
