import json
import numpy as np
import itertools
import matplotlib
from matplotlib import mlab
import matplotlib.pyplot as pp



def iterobjects(fn):
    for line in open(fn, 'r'):
        if line.startswith('#'):
            print line
            continue
        obj = json.loads(line)
        # print obj['n_states'], obj['train_lag_time']
        yield obj


def gaussianoverlap(means, vars):
    n_states, n_features = means.shape
    log_2_pi = np.log(2*np.pi)
    log_overlap = np.zeros((n_states, n_states))
    for i in range(n_states):
        for j in range(n_states):
            sigma2 = vars[i] + vars[j]
            deltamu = means[i] - means[j]
            lpr = -0.5*(log_2_pi + np.log(sigma2) + deltamu**2 / sigma2)
            log_overlap[i,j] = np.sum(lpr)

    for i in range(n_states):
        for j in range(n_states):
            log_overlap[i, j] -= 0.5*(log_overlap[i, i] + log_overlap[j, j])

    return log_overlap


for obj in iterobjects('hmm-fit-32.superpose.jsonlines'):
    transmat = np.array(obj['transmat'])
    pp.figure()
    #pp.subplot(1,3,1)
    #pp.imshow(transmat, interpolation='none')
    #pp.colorbar()

    #pp.subplot(1,3,2)
    overlap = gaussianoverlap(np.array(obj['means']), np.array(obj['vars']))
    #pp.imshow(overlap, interpolation='none')

    #pp.subplot(1,3,3)
    pp.ylabel('overlap'); pp.xlabel('transmat')
    pp.scatter(np.log(transmat).flatten(), overlap.flatten())

    #pp.colorbar()
    pp.show()
    #exit()
