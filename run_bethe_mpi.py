#!/usr/bin/env python3
import os
os.environ['OMP_NUM_THREADS'] = '1'
from mpi4py import MPI
import numpy as np
from gfmod import *

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

niwn = 60
ntau = 121
beta = 20
U = 3.0
mu = U / 2
D = 1
nwarm = 100
nmc_per_meas = 10
nmeas = 2000
nmeas_per_cpu = nmeas//comm.Get_size()
seed = 0
seedDist = 2*rank*ntau*(nmeas*nmc_per_meas+nwarm)
nloop = 30
dG = 1e-3

def print_log_msg(logtype, msg, comm):
    if comm.Get_rank() == 0:
        print(f'[{logtype}] : {msg}')

def allreduce(target, comm):
    dtype = target.dtype
    if dtype == 'float64':
        mpi_dtype = MPI.DOUBLE
    else:
        raise Exception("Check target.dtype : ", dtype)
    target_all = np.zeros_like(target)
    comm.Reduce([target, mpi_dtype], [target_all, mpi_dtype], op=MPI.SUM, root=0)
    comm.Bcast(target_all, root=0)
    return np.array(target_all/comm.Get_size(), dtype = dtype)


iwn = np.array([1j*(2*n+1)*np.pi/beta for n in range(niwn)])
Giwn = gftools.SemiCircularGreen(niwn, beta, 0, 1, D) # initial guess
tau = np.linspace(0, beta, ntau)
Gtau = gftools.GreenInverseFourier(Giwn, ntau, beta)

for niter in range(nloop):
    print_log_msg("INFO", f"iterations : {niter+1}", comm)
    giwn0 = 1/(iwn + mu - (D/2)**2 * Giwn - U/2) # In Hirsh-Fye algorithm, chemical potential of the bath Green's function is shifted by U/2.
    gtau0 = gftools.GreenInverseFourier(giwn0, ntau, beta)
    solver = qmc.hfqmc(gtau0, beta, U, seed, seedDist)
    print_log_msg("INFO", f"warm-up in MCMC ({nwarm} monte carlo steps)", comm)
    solver.do_montecarlo_step(nwarm) # warm up
    print_log_msg("INFO", f"measuring full Green's function ({nmeas_per_cpu} measurements per cpu)", comm)
    Gtau_up, err_up, Gtau_dw, err_dw = solver.meas(nmc_per_meas, nmeas_per_cpu)
    comm.Barrier()
    Gtau_new = allreduce((Gtau_up + Gtau_dw) / 2, comm) # A system is in the paramagnetic phase.
    Gtau_err = np.sqrt(allreduce(((err_up + err_dw) / 2)**2, comm)/comm.Get_size()**3)
    if np.mean(Gtau_err) > dG:
        print_log_msg("WARNING", f"Measurement error is larger than 'dG({dG})'. One should increase the # of measurements.", comm)
    dist = np.mean(np.abs(Gtau_new - Gtau))
    Gtau = np.array(Gtau_new)
    Giwn = gftools.GreenFourier(Gtau, niwn, beta)
    print_log_msg("INFO", f"MEAN(|G_new - G_old|) : {dist:.3E},  ERR(G_new) : {np.mean(Gtau_err):.3E}\n", comm)
    if dist < dG:
        print_log_msg("INFO", "CONVERGE!", comm)
        break

if rank == 0:
    np.savetxt(f"Gtau-{U}.dat", np.array([tau, Gtau, Gtau_err]).T)
    self_energy = (iwn + mu - (D/2)**2 * Giwn) - 1/Giwn # Dyson equation : \Sigma(iwn) = G0^{-1}(iwn) - G^{-1}(iwn)
    np.savetxt(f"self_energy-{U}.dat", np.array([iwn.imag, self_energy.real, self_energy.imag]).T)
