import numpy as np
import sys
sys.path.append('../GOCAM_Project/dev')

import rpy2
from rpy2.robjects.packages import importr
BiasedUrn = importr('BiasedUrn')

import utils


def get_M_wM():
    """ returns M, the number of entities in the background, and w_M, the mean size of entities in the background"""
    setID2members = utils.csv2dict('../data/setID2members.csv')
    l = []
    for s,m in setID2members.items():
        l.append(len(m))
    l = np.array(l)
    l = np.sort(l)
    num_empty_sets = np.sum(l==0)
    
    l = l[l!=0]
    mean = np.mean(l)#l[4:-4]) 1% trimmed mean?
    num_sets = len(l)
    bg = len(utils.csv2dict('../data/ID2gocam_mouse.csv'))
    M = bg-num_empty_sets
    
    w_M = np.round(((M-num_sets)+num_sets*mean)/M,decimals=2)
    return M, w_M

def make_initial_vectors(gocam2ID,setID2members, gc, M,w_M):
    """initializes counts vector (m) and weights vector (w), where each entity gets its own element in the arrays
- values in m only take on 0 (if there is no solo proteins) or 1
- values in w correspond to the weight of each element in m (weighted by the # genes in a set or 1 for solo proteins)"""
    w_gc = [1] #initialize with 1 as the weight of single proteins (irrespective of whether there are any)
    m_gc = [0] #initialize with 0 single proteins
    num_protein = 0
    for i in gocam2ID.get(gc):
        if "sset:" in i:
            w_i = len(setID2members.get(i))
            w_gc.append(w_i)
            m_gc.append(1)
        else:
            num_protein+=1
    m_gc[0] = num_protein
    m_gc.append(M-np.sum(m_gc)) #entities not in the gocam (roughly)
    w_gc.append(w_M) #weight for entities not in the gocam (all weighted as w_M (the mean))
    return w_gc, m_gc


def make_new_vectors(w_gc,m_gc,M,w_M):
    """compress the m and w vectors by grouping elements according to their weights
- w is the ordered set of unique weights for entities of the gocam + the background bin
- m[i] is the number of entities in the pathway with the weight specified in w[i] + the background bin"""
    w_temp = w_gc[:-1]
    if w_temp[0] != 1:
        print('Possible bug: w_temp[0] != 1',w_temp)
        
    w_new, m_temp = np.unique(w_temp, return_counts=True)
    m_temp[0]=m_gc[0] #w_gc and m_gc have weight 1 as w_gc[0] and the number of single proteins as m_gc[0]
    m_new = np.append(m_temp,np.array([M-np.sum(m_temp)]))
    w_new = np.append(np.unique(w_temp),np.array([w_M]))
    return w_new, m_new




def ncHGT_sf(XT,m,N,w):
    """survival function, sums PMF for all possibilities where K >= k by calling BiasedUrn"""
    #l = len(XT)/len(m)
    if len(XT) == 0:
        print('len(XT) = 0')
        return -1
    pval = 0
    #np.seterr(under='warn')
    for i in range(len(XT)):
        x = rpy2.robjects.IntVector(XT[i])
        pval = pval + BiasedUrn.dMFNCHypergeo(x,m,N,w, precision = 1e-10)[0]
    return pval


def enumerate_possibilities(m_new,i,prev_array):
    """enumerate all possible counts vectors"""
    first = True
    for j in range(m_new[i]+1):
        xt = prev_array.copy()
        xt[0][i] = j
        
        #recursion
        if (i < len(m_new)-1):
            xt = enumerate_possibilities(m_new, i+1, xt) #will return matrix (array of arrays)
            
        #combining results into matrix
        if not first:
            XT = np.concatenate([XT,xt], axis = 0)
        else:
            XT = xt
            first = False
    return XT


def do_ncHGT(k,gc,M,N):
    setID2members = utils.csv2dict('../data/setID2members.csv')
    gocam2ID = utils.csv2dict('../data/gocam2ID_mouse.csv')
    
    M, w_M = get_M_wM()
    
    #make weight (w) and bin size (m) vectors where each entity in the gocam gets its own entry
    w_in, m_in = make_initial_vectors(gocam2ID, setID2members, gc, M,w_M)
    
    #update m and w vectors by grouping sets of the same size
    w_new , m_new= make_new_vectors(w_in,m_in,M,w_M)

    #make XT matrix, an enumeration of all possible arangements of balls in bins based on m_new and w_new
    m_gc = m_new[:-1] #don't pass the background bin to XT
    XT = enumerate_possibilities(m_gc,0,np.zeros(shape=(1,len(m_gc))))
    
    #filter XT to only include the region of the sample space >= k (which is what we want to sum probabilities over)
    mask1 = (np.sum(XT, axis=1) >= k)
    XT = XT[mask1]

    #filter XT to ensure that more than N entities are not picked
    mask2 = (np.sum(XT, axis=1) <= N)
    XT = XT[mask2]

    #add the remaining entities to the m+1th bin (non gocam bin)
    x_mp1_vec = N- np.sum(XT, axis = 1) #number of balls to be drawn from the last bin (the non-gocam background)
    XT = np.concatenate((XT,x_mp1_vec.reshape(len(x_mp1_vec),1)), axis = 1)

    m = rpy2.robjects.IntVector(m_new)
    w = rpy2.robjects.FloatVector(w_new)
    pval = ncHGT_sf(XT,m,N,w)
    return pval


