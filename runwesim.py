import numpy as np
import os

def program_paths():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    cwesis_path = dir_path +'/cwesis.exe'
    netcraft_path = dir_path +'/netcraft.py'
    slurm_path = dir_path + '/slurm.serjob'
    return cwesis_path,netcraft_path,slurm_path


def par_paths(net_num):
    data_path = os.getcwd() + '/'
    path_adj_in = data_path + 'Adjin_{}.txt'.format(net_num)
    path_adj_out = data_path + 'Adjout_{}.txt'.format(net_num)
    path_parameters = data_path + 'cparameters_{}.txt'.format(net_num)
    return path_adj_in,path_adj_out,path_parameters


def run_program(N,lam,number_of_networks,k,x,Num_inf,Alpha,Beta_avg,eps_din,eps_dout,network_type,sims,relaxation_time,
                it,tau,jump,new_trajcetory_bin,Istar):
    cwesis_path, netcraft_path, slurm_path = program_paths()
    path_all_results = os.getcwd()
    for n,l,net,degree,inf,gamma,b,epsin,epsout,net_type,s,rxt,runs,t,j,ninf,cutoff,avg_inf in zip(N,lam,
               number_of_networks,k,x,Num_inf,Alpha,Beta_avg,eps_din,eps_dout,network_type,sims,relaxation_time,it,
                tau,jump,new_trajcetory_bin,Istar):
        os.chdir(path_all_results)
        foldername = 'prog_{}_N{}_k_{}_R_{}_tau_{}_it_{}_jump_{}_new_trajcetory_bin_{}_sims_{}_epsin_{}_epsout_{}'.format(
            net_type, n, degree, l, t, runs, j, cutoff, s, epsin, epsout)
        os.mkdir(foldername)
        os.chdir(foldername)
        data_path = os.getcwd()
        for i in range(net):
            path_adj_in,path_adj_out,path_parameters = par_paths(i)
            parameters = np.array([n,s,runs,degree,inf,l,j,ninf,gamma,b,net,t,cutoff,net_type,epsin,epsout,avg_inf])
            np.save('parameters_{}.npy'.format(i), parameters)
            os.system('{} {} {} {} {}'.format(slurm_path,netcraft_path,cwesis_path,path_adj_in,path_adj_out,path_parameters))
            os.chdir(data_path)

if __name__ == '__main__':
    # Network parameters
    N = np.array([100,200]) # number of nodes
    lam = np.array([1.6,1.6]) # The reproduction number
    number_of_networks = np.array([10,10])
    k = np.array([50,50]) # Average number of neighbors for each node
    eps_din,eps_dout = np.array([0.0,0.0]),np.array([0.0,0.0]) # The normalized std (second moment divided by the first) of the network
    x = 0.2*np.ones(len(N)) # intial infection percentage
    Num_inf = x*N # Number of initially infected nodes
    Alpha = 1.0*np.ones(len(N)) # Recovery rate
    Beta_avg = Alpha * lam / k # Infection rate for each node
    network_type = np.repeat('h',len(N))

    # WE simulation parameters
    sims = 200*np.ones(len(N)) # Number of simulations at each step
    relaxation_time  = 3*np.ones(len(N))
    it = 70*np.ones(len(N)) # Number of iterations WE program runs
    tau = 1.0*np.ones(len(N))
    jump = 15*np.ones(len(N))
    new_trajcetory_bin = 15*np.ones(len(N))


    # Number of infected fixed point
    # y1star=(-2*eps_din*(1 + eps_dout*eps_din)+ lam*(-1 + eps_din)*(1 + (-1 + 2*eps_dout)*eps_din)+ np.sqrt(lam**2 +eps_din*(4*eps_din +lam**2*eps_din*(-2 +eps_din**2) +4*eps_dout*(lam -(-2 + lam)*eps_din**2) +4*eps_dout**2*eps_din*(lam -(-1 + lam)*eps_din**2))))/(4*lam*(-1 +eps_dout)*(-1 +eps_din)*eps_din)
    # y2star=(lam + eps_din*(-2 + 2*lam +lam*eps_din+ 2*eps_dout*(lam +(-1 + lam)*eps_din)) -np.sqrt(lam**2 +eps_din*(4*eps_din +lam**2*eps_din*(-2 +eps_din**2) +4*eps_dout*(lam -(-2 + lam)*eps_din**2) +4*eps_dout**2*eps_din*(lam -(-1 + lam)*eps_din**2))))/(4*lam*(1 +eps_dout)*eps_din*(1 + eps_din))
    # Istar = (y1star +y2star)*N
    Istar = (1 - 1/lam) * N

    # Run the slurm program
    run_program(N,lam,number_of_networks,k,x,Num_inf,Alpha,Beta_avg,eps_din,eps_dout,network_type,sims,relaxation_time,it,
                tau,jump,new_trajcetory_bin,Istar)

