#!/usr/bin/env python3
import os
os.environ['OMP_NUM_THREADS'] = '1'

import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray
from mpi4py import MPI
from loguru import logger

from gfmod import gftools, qmc


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

beta = 20  # inverse temperature
U = 3.0  # on-site interaction
V = 1.0  # hybridization strength
D = 1.0  # half-bandwidth of semi-circular DOS
mu = -3.0  # chemical potential
niwn = 120  # cut-off number of Matsubara frequency
ntau = int(4 * beta * U) + 1  # Imaginary time is discritized in 'ntau' slice.
nmc_per_meas = 10  # number of Monte-Carlo steps for each measurement
nwarm = 100  # number of Monte-Carlo steps for warm-up
nmeas = 5000  # number of measurements
nmeas_per_cpu = nmeas // comm.Get_size()
seed = 0  # random seed
seedDist = 2 * rank * ntau * (nmeas * nmc_per_meas + nwarm)  # Each processor has different random seed.
tau = np.linspace(0, beta, ntau)  # imaginary time (0~beta)


def print_log_msg(msg: str, comm) -> None:
    if comm.Get_rank() == 0:
        logger.info(f'{msg}')
    comm.Barrier()


def allreduce(target: str, comm) -> NDArray[np.float64]:
    dtype = target.dtype
    if dtype == 'float64':
        mpi_dtype = MPI.DOUBLE
    else:
        raise Exception("Check target.dtype : ", dtype)
    target_all = np.zeros_like(target)
    comm.Reduce([target, mpi_dtype], [target_all, mpi_dtype], op=MPI.SUM, root=0)
    comm.Bcast(target_all, root=0)
    return np.array(target_all / comm.Get_size(), dtype=dtype)


ax = plt.axes((0, 0, 0.4, 0.4))
ax.set_xlabel(r'$\beta$', fontsize=10)
ax.set_ylabel(r'$-G(\tau)$', fontsize=10)
ax.set_xlim(0, beta)
ax.set_ylim(0, 1)
ax.tick_params(which='both', direction='in')

giwn = gftools.SemiCircularGreen(niwn, beta, mu, V, D)
gtau = gftools.GreenInverseFourier(giwn, ntau, beta, (1, 0, 0,))
solver = qmc.hfqmc(gtau, beta, U, seed, seedDist)
print_log_msg(f"warm-up in MCMC ({nwarm} monte carlo steps)", comm)
solver.do_montecarlo_step(nwarm)  # warm up
print_log_msg(f"measuring full Green's function ({nmeas_per_cpu} measurements per cpu)", comm)
Gtau_up, err_up, Gtau_dw, err_dw = solver.meas(nmc_per_meas, nmeas_per_cpu)
comm.Barrier()
Gtau = allreduce((Gtau_up + Gtau_dw) / 2, comm)  # due to Z2 symmetry
Gtau_err = np.sqrt(allreduce(((err_up + err_dw) / 2)**2, comm) / comm.Get_size()**3)
if comm.Get_rank() == 0:
    ax.errorbar(tau, -Gtau, Gtau_err,
        linestyle='-',
        label=r'$\mu={}$'.format(mu)
    )
comm.Barrier()

if comm.Get_rank() == 0:
    ax.legend(loc='best', fontsize=10)
    plt.savefig('anderson_model.png', dpi=600, bbox_inches='tight')
