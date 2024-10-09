# This program run multiple WE simulation on multiple networks with different parameters
import os
import numpy as np

if __name__ == '__main__':

    # Netwrok parameters

    # N = [300,400,500,600,700,800,900,1000,1100,1200,1300,1400]
    N = 5000
    prog = 'bd'
    lam = 1.2
    # lam = 1+np.logspace(-2,0,9)
    # lam = np.array([1.5,1.6,1.7,1.8])
    eps_din = 0.5
    eps_dout = 0.5
    #eps_din = [0.01,0.04,0.06,0.08,0.1,0.14,0.18,0.2,0.25,0.3,0.4,0.5,0.6]
    #eps_dout = [0.01,0.04,0.06,0.08,0.1,0.14,0.18,0.2,0.25,0.3,0.4,0.5,0.6]
    skewness = [-0.6,-0.5,-0.4,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6]
    # correlation = 0.3
    number_of_networks = 10
    # k = [50]
    k= 50

    # We simulation parameters

    sims = 500
    tau = 1.0
    # tau = np.linspace(0.1,2.0,20)
    it = 70
    jump = 1
    new_trajectory_bin = 2
    error_graphs = False

    # Parameters that don't change

    relaxation_time = 20
    x = 0.2
    Alpha = 1.0

    # Paths needed to run the program
    dir_path = os.path.dirname(os.path.realpath(__file__))
    slurm_path = dir_path +'/slurm.serjob python3'
    program_path = dir_path +'/runwesim.py'
    loop_over = skewness

    for i in loop_over:
        error_graphs_flag = '--error_graphs' if error_graphs else ''

        command = (f'{slurm_path} {program_path} --N {N} --prog {prog} --lam {lam} --eps_din {eps_din} '
                   f'--eps_dout {eps_dout} --skewness {i} --number_of_networks {number_of_networks} '
                   f'--k {k} {error_graphs_flag} --sims {sims} --tau {tau} --it {it} --jump {jump} '
                   f'--new_trajectory_bin {new_trajectory_bin} --relaxation_time {relaxation_time} --x {x} '
                   f'--Alpha {Alpha}')
        os.system(command)