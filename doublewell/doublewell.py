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


def exact_solution(n_grid, plot_eigfunctions=True):
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
        pp.savefig('../figures/doublewell-eigenfunctions-msm.png')
        exit(1)


    eigs = np.sort(np.real(np.linalg.eigvals(transmat)))
    timescales =  -1 / np.log(eigs)
    return timescales


def msm_solution(x, n_grid, lag_time):
    x = x[::lag_time]

    countsmat = np.zeros((n_grid, n_grid), dtype=np.float)
    transmat = np.zeros((n_grid, n_grid), dtype=np.float)    
    grid = np.linspace(-np.pi, np.pi, n_grid)

    labels = np.zeros_like(x)
    for i in range(len(x)):
        labels[i] = np.argmin((x[i] - grid)**2)
    for i in range(len(x)-1):
        countsmat[labels[i], labels[i+1]] += 1
    countsmat += 0.001

    for i in range(n_grid):
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

def msm_timescale(trajectories, lag_times):
    all_timescales = np.zeros((len(trajectories), len(lag_times)))
    for i, x in enumerate(trajectories):
        all_timescales[i] = [msm_solution(x, 2, lag_time)[-2] for lag_time in lag_times]
    return  np.mean(all_timescales, axis=0), np.std(all_timescales, axis=0) / np.sqrt(len(trajectories))

def ghmm_timescale(trajectories, lag_times):
    all_timescales = np.zeros((len(trajectories), len(lag_times)))
    for i, x in enumerate(trajectories):
        all_timescales[i] = [GaussianFusionHMM(2, n_features=1, fusion_prior=0).fit([x[::l].reshape(-1,1)]).timescales_()[-1]*l for l in lag_times]
    return  np.mean(all_timescales, axis=0), np.std(all_timescales, axis=0) / np.sqrt(len(trajectories))


exact_solution(300)

lag_times = [25, 50, 100, 200, 250, 500, 750, 1000, 1300, 1600, 1950]
msm_mean, msm_std = msm_timescale(trajectories, lag_times)
hmm_mean, hmm_std = ghmm_timescale(trajectories, lag_times)

pp.ion()
pp.figure(figsize=(8,10))
pp.subplot(313)
pp.ylabel('Relaxation Time')
pp.xlabel('Lag Time')
pp.errorbar(lag_times, msm_mean, msm_std, ls='--', label='2-state MSM', color=colors[0], lw=2)
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
ip.embed()