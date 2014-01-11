import os
import IPython as ip
import numpy as np
from mdtraj import io
from mixtape.ghmm import GaussianFusionHMM
import matplotlib.pyplot as pp
colors = ['#348ABD', '#7A68A6', '#A60628', '#467821', '#CF4457', '#188487', '#E24A33']

DIFFUSION_CONST = 1
DT = 1e-3
SQRT_2D_DT = np.sqrt(2 * DIFFUSION_CONST * DT)
grad_potential = lambda x : -2 * np.sin(2 * x)


def reflect_pbc(x, min, max):
    if x > max:
        return 2*max - x
    if x < min:
        return 2*min - x
    return x


def propagate(n_steps):
    import time; start = time.time()
    n_steps = int(n_steps)

    max = np.pi
    min = -np.pi
    rand = np.random.randn(n_steps)
    x = np.zeros(n_steps+1)
    for i in range(n_steps):
        x_i_plus_1 = x[i] -DT * grad_potential(x[i]) + SQRT_2D_DT*rand[i]
        x[i+1] = reflect_pbc(x_i_plus_1, -np.pi, np.pi)
        # reflecting PBCs

    print '%d steps/s' % (n_steps / (time.time() - start))
    return x


def exact_solution(n_grid, plot_eigfunctions=False):
    ONE_OVER_SQRT_2PI = 1.0 / (np.sqrt(2*np.pi))
    normalpdf = lambda x : ONE_OVER_SQRT_2PI * np.exp(-0.5 * (x*x))

    grid = np.linspace(-np.pi, np.pi, n_grid)
    width = grid[1]-grid[0]
    transmat = np.zeros((n_grid, n_grid))
    for i, x_i in enumerate(grid):
        for offset in range(-(n_grid-1), n_grid):
            x_j = x_i + (offset * width)
            j = reflect_pbc(i+offset, 0, n_grid-1)

            # What is the probability of going from x_i to x_j in one step?
            diff = (x_j - x_i + DT * grad_potential(x_i)) / SQRT_2D_DT
            transmat[i, j] += normalpdf(diff)
        #print transmat[i, :]
        transmat[i, :] =  transmat[i, :] / np.sum(transmat[i, :])


    eigvalues, eigvectors = np.linalg.eig(transmat)
    eigsort = np.argsort(np.real(eigvalues))[::-1]
    eigvectors = eigvectors[:, eigsort]
    eigvalues = eigvalues[eigsort]
    
    if plot_eigfunctions:
        pp.title('Double well transfer operator 1st eigenfunction')
        eigvectors[:, 1] /= np.max(eigvectors[:, 1])
        
        pp.plot(grid, eigvectors[:, 1], label='Exact Soln.', lw=3)
        pp.plot([-np.pi, 0, 0, np.pi], [0.85, 0.85, -0.85, -0.85], label='2-state MSM', lw=2)

        xx = np.linspace(-np.pi, np.pi, 10+1)
        vv = np.zeros(10+1)
        for i in range(10):
            center = (xx[i] + xx[i+1]) / 2
            vv[i] = eigvectors[np.argmin((grid-center)**2), 1]
        for i in range(1, 10):
            pp.plot([xx[i], xx[i]], [vv[i-1], vv[i]], c='k', lw=2)
        for i in range(10):
            pp.plot([xx[i], xx[i+1]], [vv[i], vv[i]], c='k', lw=2)
        pp.plot([0,0], [0,0], c='k', lw=2, label='10-state MSM')
        #pp.yticks([])
        pp.legend()
        pp.ylim(-1.2, 1.2)
        pp.xlim(-np.pi, np.pi)
        print 'saving eigenfunctions'
        pp.savefig('doublewell-eigenfunctions-msm.png')
        exit(1)


    eigs = np.sort(np.real(np.linalg.eigvals(transmat)))
    timescales =  -1 / np.log(eigs)
    return timescales


def msm_solution(x, k, lag_time, discretization='grid'):
    x = x[::lag_time]

    countsmat = np.zeros((k, k), dtype=np.float)
    transmat = np.zeros((k, k), dtype=np.float)    

    if discretization == 'grid':
        grid = np.linspace(-np.pi, np.pi, k)
        labels = np.zeros_like(x)
        for i in range(len(x)):
            labels[i] = np.argmin((x[i] - grid)**2)
    elif hasattr(discretization, '__call__'):
        # manual discretization
        labels = discretization(x)
    else:
        raise NotImplementedError()

    for i in range(len(x)-1):
        countsmat[labels[i], labels[i+1]] += 1
    countsmat += 0.001

    for i in range(k):
        transmat[i, :] =  countsmat[i, :] / np.sum(countsmat[i, :])

    eigs = np.sort(np.real(np.linalg.eigvals(transmat)))
    timescales =  -lag_time / np.log(eigs)
    return timescales

trajectories = []
for i in range(10):
    fn = 'trajectory-%d.h5' % i
    if os.path.exists(fn):
        print 'loading %s' % fn
        trajectories.append(io.loadh(fn)['arr_0'])
    else:
        x = propagate(5e5)
        io.saveh(fn, x)
        print 'saving %s' % fn
        trajectories.append(x)

def msm_timescale(trajectories, lag_times, n_states=2, discretization='grid'):
    all_timescales = np.zeros((len(trajectories), len(lag_times)))
    for i, x in enumerate(trajectories):
        all_timescales[i] = [msm_solution(x, n_states, lag_time, discretization=discretization)[-2] for lag_time in lag_times]
    return  np.mean(all_timescales, axis=0), np.std(all_timescales, axis=0) / np.sqrt(len(trajectories))

def ghmm_timescale(trajectories, lag_times, n_states=2):
    all_timescales = np.zeros((len(trajectories), len(lag_times)))
    for i, x in enumerate(trajectories):
        all_timescales[i] = [GaussianFusionHMM(n_states, n_features=1, fusion_prior=0).fit([x[::l].reshape(-1,1)]).timescales_()[-1]*l for l in lag_times]
    return  np.mean(all_timescales, axis=0), np.std(all_timescales, axis=0) / np.sqrt(len(trajectories))


def plot_stepfunction(edges, values, c='k', lw=2, label=None):
    assert len(edges)-1 == len(values)
    for i in range(len(edges)-1):
        pp.plot([edges[i], edges[i+1]], [values[i], values[i]], color=c, lw=lw)
    for i in range(len(values)-1):
        pp.plot([edges[i+1], edges[i+1]], [values[i], values[i+1]], color=c, lw=lw)
    if label is not None:
        pp.plot([0], [0], lw=lw, c=c, label=label)
    pass

def main1():
    colors = ['r', 'b']
    from scipy.cluster.vq import kmeans, vq
    #centroids, distortion = kmeans(trajectories[0], 4)
    centroids = np.sort(np.array([-1.11630864, -2.14613473,  2.1049366,   1.10232638]))
    print 'centroids', centroids

    manual_bin_edges = [-1, 0, 1]
    manual = lambda x: np.digitize(x, manual_bin_edges)
    means = lambda x: vq(x, centroids)[0]
    
    
    lag_times = [25, 50, 100, 200, 250, 500, 750, 1000, 1300, 1600, 1950]
    #lag_times = [1000,]
    timescale_kmeans_avg, timescale_kmeans_std = msm_timescale(trajectories, lag_times, n_states=4, discretization=means)
    timescale_manual_avg, timescale_manual_std = msm_timescale(trajectories, lag_times, n_states=4, discretization=manual)

    pp.figure(figsize=(10,8))
    pp.subplot(222)
    pp.errorbar(lag_times, timescale_kmeans_avg, timescale_kmeans_std,
                ls='-', label='kmeans "optimal" states', color=colors[0], lw=2)
    pp.errorbar(lag_times, timescale_manual_avg, timescale_manual_std,
                ls='-', label='manual states', color=colors[1], lw=2)
    pp.ylabel('Relaxation Timescales')
    pp.xlabel('Lag Time')
    exact = exact_solution(300)[-2]
    pp.plot([0, max(lag_times)], [exact, exact], lw=2, c='k', label='Exact')
    pp.legend(loc=4)#, prop={'size': 14})
    
    
    pp.subplot(221)
    pp.hist(np.array(trajectories).reshape(-1), bins=50, log=True, color='grey', label='Samples')
    pp.ylabel('Frequency')
    pp.legend(loc=2)
    pp.twinx()
    pp.plot(np.linspace(-np.pi, np.pi), 1+np.cos(2*np.linspace(-np.pi, np.pi)), c='k', label='Potential', lw=2)
    pp.xlabel('Position')
    pp.ylabel('E / kT')
    pp.xlim(-np.pi, np.pi)
    pp.legend(loc=1)

    
    
    pp.subplot(223)
    pp.plot(np.linspace(-np.pi, np.pi), 1+np.cos(2*np.linspace(-np.pi, np.pi)), '-', c='k', lw=1)
    for e in manual_bin_edges:
        pp.plot([e,e], [0, 2], '--', lw=3, c=colors[1])
    pp.plot([0], [0], lw=2, c=colors[1], label='manual states')
    pp.text(-2, 1, '1', color=colors[1], size=32, ha='center')
    pp.text(-0.5, 1, '2', color=colors[1], size=32, ha='center')
    pp.text(0.5, 1, '3', color=colors[1], size=32, ha='center')
    pp.text(2, 1, '4', color=colors[1], size=32, ha='center')
    pp.yticks([])
    pp.xlabel('Position')
    pp.legend(loc=4)#, prop={'size': 14})
    pp.xlim(-np.pi, np.pi)

    pp.subplot(224)
    edges_from_centroids = lambda l: [0.5*(l[i] + l[i+1]) for i in range(len(l)-1)]
    pp.plot(np.linspace(-np.pi, np.pi), 1+np.cos(2*np.linspace(-np.pi, np.pi)), '-', c='k', lw=1)
    for e in edges_from_centroids(centroids):
        pp.plot([e,e], [0, 2], '--', lw=3, c=colors[0])
    pp.plot([0], [0], lw=2, c=colors[0], label='kmeans states')
    for i in range(len(centroids)):
        pp.text(centroids[i], 1, str(i+1), color=colors[0], size=32, ha='center')
    pp.yticks([])
    pp.xlabel('Position')
    pp.legend(loc=4)#, prop={'size': 14})
    pp.xlim(-np.pi, np.pi)

    #pp.tight_layout()
    #pp.subplots_adjust(wspace=0.3)
    ip.embed()
    exit(1)

if __name__ == '__main__':
    pp.ion()
    main1()

#exact_solution(300)

lag_times = [25, 50, 100, 200, 250, 500, 750, 1000, 1300, 1600, 1950]
msm_mean, msm_std = msm_timescale(trajectories, lag_times, n_states=4, discretization='kmeans')
msm_mean2, msm_std2 = msm_timescale(trajectories, lag_times, n_states=4, discretization=lambda x: np.digitize(x, [-1, 0, 1]))

hmm_mean, hmm_std = ghmm_timescale(trajectories, lag_times)

pp.ion()
pp.figure(figsize=(8,10))
pp.subplot(313)
pp.ylabel('Relaxation Time')
pp.xlabel('Lag Time')
pp.errorbar(lag_times, msm_mean, msm_std, ls='--', label='4-state MSM kmeans', color=colors[0], lw=2)
pp.errorbar(lag_times, msm_mean2, msm_std2, ls='--', label='4-state MSM manual', color=colors[1], lw=2)

#pp.fill_between(lag_times, msm_mean-msm_std, msm_mean+msm_std, color=colors[0], alpha=0.2)

#pp.errorbar(lag_times, hmm_mean, hmm_std, ls='-', label='2-state gHMM', color=colors[1], lw=2)
#pp.fill_between(lag_times, hmm_mean-hmm_std, hmm_mean+hmm_std, color=colors[1], alpha=0.2)

exact = exact_solution(300)[-2]
pp.plot([0, max(lag_times)], [exact, exact], lw=2, c='k', label='Exact')
pp.legend(loc=4, prop={'size': 14})
pp.xlim(0, 2000)
pp.ylim(0, 8000)


pp.subplot(311)
pp.hist(np.array(trajectories).reshape(-1), bins=50, log=True, label='Samples')
pp.ylabel('Frequency')
pp.legend(loc=2)
pp.twinx()
pp.plot(np.linspace(-np.pi, np.pi), 1+np.cos(2*np.linspace(-np.pi, np.pi)), c='k', label='Potential', lw=2)
pp.xlabel('Position')
pp.ylabel('E / kT')
pp.xlim(-np.pi, np.pi)
pp.legend(loc=1)


pp.subplot(312)
pp.plot(np.arange(len(trajectories[0])), trajectories[0], c='k')
pp.xlabel('Time [steps]')
pp.ylabel('Position')

#pp.show()
print '#'*80
ip.embed()
