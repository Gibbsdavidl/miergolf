__author__ = 'David L Gibbs'

import scipy.sparse as sp
import scipy.sparse.linalg as lin
import numpy as np

def weightsum(nodes,tup):
    totwt = 0.0
    for ti in tup:
        totwt += nodes[ti][2]
    return(totwt)


def scoreSoln(soln, s, smat, nodes):
    # soln -- the solutions set S
    # smat  -- nxn sparse matrix
    # s     -- the program state
    wt = weightsum(nodes,soln)
    n = (smat.shape[0])
    ts = [i for i in xrange(n) if i not in soln] # set T
    idn = sp.eye(len(ts),len(ts))
    pst = subMatrix(soln, ts, smat) # prob of S to T
    pts = subMatrix(ts, soln, smat) # prob of T to S
    ptt = subMatrix(ts, ts, smat)    # prob of T to T
    lap = sp.csc_matrix(idn-ptt)
    pst_t = sp.csc_matrix(pst.transpose())
    lap_t = sp.csc_matrix(lap.transpose())
    if s["mode"] == "both":
        return(scoreBoth(s,lap,pts,lap_t,pst_t,wt))
    elif s["mode"] == "tx":
        return(scoreTX(s,lap_t,pst_t,wt))
    elif s["mode"] == "rx":
        return(scoreRX(s,lap,pts,wt))
    else:
        print "ScoreSoln Error! mode must be rx, tx, or both."
        sys.exit(1)


def scoreBoth(s,lap,pts,lap_t,pst_t,wt):
    try:
        f = lin.spsolve(lap, pts)
        h = lin.spsolve(lap_t, pst_t)
    except:
        # a singular matrix ... must solve each vector separately
        vecs = pts.shape[1]
        vlen = pts.shape[0]
        #print(vecs)
        fsolns = []
        hsolns = []
        for i in range(vecs):
            pts2   = pts[:,i].todense()
            pst_t2 = pst_t[:,i].todense()
            f1 = lin.bicgstab(lap, pts2)[0]
            h1 = lin.bicgstab(lap_t, pst_t2)[0]
            fsolns.append(f1)
            hsolns.append(h1)
        f = np.matrix(fsolns)
        h = np.matrix(hsolns)
        f = f.transpose()
        h = h.transpose()
    if type(f) == type(np.array([])): # came back as an array
        fh = f+h
        score = fh.sum()
        touch = (fh > s["tx"]).sum()
    else: # came back as a sparse matrix
        fsum = np.array(f.sum(1)).flatten()
        hsum = np.array(h.sum(1)).flatten()
        fh = fsum + hsum
        touch = sum(fh > s["tx"])
        #  best score; best touch #
    return((touch+wt, touch))



def scoreRX(s,lap,pts,wt):
    f = lin.spsolve(lap, pts)
    if type(f) == type(np.array([])): # came back as an array
        ftouch = sum(f > s["rx"])
    else: # came back as a sparse matrix
        fsum = np.array(f.sum(1)).flatten()
        ftouch = sum(fsum > s["rx"])
    return((ftouch+wt, ftouch))


def scoreTX(s, lap_t, pst_t, wt):
    try:
        h = lin.spsolve(lap_t, pst_t)
    except:
        # a singular matrix ... must solve each vector separately
        vecs = pst_t.shape[1]
        hsolns = []
        for i in range(vecs):
            pst_t2 = pst_t[:,i].todense()
            h1 = lin.bicgstab(lap_t, pst_t2)[0]
            hsolns.append(h1)
        h = np.matrix(hsolns)
        h = h.transpose()
    if type(h) == type(np.array([])): # came back as an array
        htouch = sum(h > s["tx"])
    else: # came back as a sparse matrix
        hsum = np.array(h.sum(1)).flatten()
        htouch = sum(hsum > s["tx"])
    return((htouch+wt, htouch))


def scoreMats(soln, s, smat):
    # ss    -- the solutions set S
    # smat  -- nxn sparse matrix
    # ts    -- the set T
    n = (smat.shape[0])
    ts = [i for i in xrange(n) if i not in soln]
    idn = sp.eye(len(ts))
    pst = subMatrix(soln, ts, smat)
    pts = subMatrix(ts, soln, smat)
    ptt = subMatrix(ts, ts, smat)
    lap = sp.csc_matrix(idn-ptt)
    pst_t = sp.csc_matrix(pst.transpose())
    lap_t = sp.csc_matrix(lap.transpose())
    f = lin.spsolve(lap, pts)
    h = lin.spsolve(lap_t, pst_t)
    return((f,h))


def subMatrix(rows, cols, A):
    a1 = A.tocsc()[:,cols]
    a2 = a1.tocsr()[rows,:]
    return(a2)
    #return(A.tocsr()[rows,:].tocsc()[:,cols])


def scoreMax(solns, s):
    sortedSolns = sorted(solns, key=lambda x: x[1][1], reverse=True)
    if (sortedSolns[0][1][1] > sortedSolns[1][1][1]):
        # if the top two scores are the same, then we need to consider the edge weights
        solns = [ x for x in solns if x[1][1] == sortedSolns[0][1][1] ] # get out the solns that have the same top touch
        sortedSolns = sorted(solns, key=lambda x: x[1][0], reverse=True) # sort by the score within that group.
    return(sortedSolns[0])
