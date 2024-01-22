import pickle
import numpy as np
import os
# import we_sis_network_extinction
import rand_networks
import networkx as nx
import csv
import pickle

def export_parameters_to_csv(parameters,network_number):
    name_parameters = 'cparameters_{}.txt'.format(network_number)
    # N, sims, it, k, x, lam, jump, Alpha,Beta,number_of_networks, tau, mf_solution ,eps_din, eps_dout, new_trajcetory_bin, prog, Beta_avg,dir_path = parameters
    # cparameters=[N, sims, it, k, x, lam, jump, Alpha,Beta,number_of_networks, tau, mf_solution ,new_trajcetory_bin, prog, Beta_avg,dir_path]
    f =open(name_parameters,'+a')
    with f:
        writer = csv.writer(f)
        writer.writerow(parameters)
    f.close()

def export_network_to_csv(G,netname):
    # Open a CSV file for writing incoming neighbors
    with open('Adjin_{}.txt'.format(netname), 'w', newline='') as incoming_file:
        # Create a CSV writer
        incoming_writer = csv.writer(incoming_file)
        # Iterate over all nodes in the graph
        for node in np.sort(G):
            # Get the incoming neighbors of the current node
            incoming_neighbors = list(G.predecessors(node))
            # Get the degree of the current node
            degree = G.in_degree[node]
            # Write a row to the CSV file for the current node
            joint = np.concatenate(([degree],incoming_neighbors),axis=0)
            incoming_writer.writerow(joint)
            # incoming_writer.writerow([degree])
            # for node in incoming_neighbors: incoming_writer.writerow([node])
    with open('Adjout_{}.txt'.format(netname), 'w', newline='') as outgoing_file:
        # Create a CSV writer
        outgoing_writer = csv.writer(outgoing_file)
        # Iterate over all nodes in the graph
        for node in np.sort(G):
            # Get the incoming neighbors of the current node
            outgoing_neighbors = list(G.successors(node))
            # Get the degree of the current node
            degree = G.out_degree[node]
            # Write a row to the CSV file for the current node
            joint = np.concatenate(([degree],outgoing_neighbors),axis=0)
            outgoing_writer.writerow(joint)

def act_as_main(foldername,parameters,Istar):
    # This program will run the we_sis_network_extinction.py on the laptop\desktop
    program_path = '/home/elad/we_sis_networks_extinction_cpp/cwesis.exe'
    os.mkdir(foldername)
    os.chdir(foldername)
    dir_path = os.path.dirname(os.path.realpath(__file__)) +'/' + foldername +'/'
    N,sims,it,k,x,lam,jump,Num_inf,Alpha,number_of_networks,tau,eps_din,eps_dout,new_trajcetory_bin,prog,Beta_avg = parameters
    if prog == 'bd':
        # G = nx.complete_graph(N)
        d1_in, d1_out, d2_in, d2_out = int(int(k) * (1 - float(eps_din))), int(int(k) * (1 - float(eps_dout))), int(int(k) * (1 + float(eps_din))), int(
            int(k) * (1 + float(eps_dout)))
        Beta = float(Beta_avg) / (1 + float(eps_din) * float(eps_dout))  # This is so networks with different std will have the reproduction number
        parameters = np.array([N,sims,it,k,x,lam,jump,Alpha,Beta,number_of_networks,tau,Istar,new_trajcetory_bin,dir_path,prog,dir_path])
        np.save('parameters.npy',parameters)
    # dir_path = os.path.dirname(os.path.realpath(__file__))
    for i in range(int(number_of_networks)):
        if prog=='bd':
            G = rand_networks.random_bimodal_directed_graph(int(d1_in), int(d1_out), int(d2_in), int(d2_out), int(N))
            parameters = np.array([N,sims,it,k,x,lam,jump,Alpha,Beta,i,tau,Istar,new_trajcetory_bin,dir_path,prog,dir_path])
        else:
            G = rand_networks.configuration_model_directed_graph(prog, float(eps_din), float(eps_dout), int(k), int(N))
            k_avg_graph = np.mean([G.in_degree(n) for n in G.nodes()])
            Beta_graph = float(lam)/k_avg_graph
            eps_in_graph = np.std([G.in_degree(n) for n in G.nodes()])/k_avg_graph
            eps_out_graph = np.std([G.out_degree(n) for n in G.nodes()])/k_avg_graph
            Beta = Beta_graph / (1 + np.sign(float(eps_din))*eps_in_graph * np.sign(float(eps_dout))* eps_out_graph)
            parameters = np.array([N,sims,it,k,x,lam,jump,Alpha,Beta,i,tau,Istar,new_trajcetory_bin,dir_path,prog,dir_path])
            np.save('parameters_{}.npy'.format(i), parameters)
        infile = 'GNull_{}.pickle'.format(i)
        with open(infile,'wb') as f:
            pickle.dump(G,f,pickle.HIGHEST_PROTOCOL)
        # nx.write_gpickle(G, infile)
        export_network_to_csv(G, i)
        export_parameters_to_csv(parameters,i)
        path_adj_in = dir_path + 'Adjin_{}.txt'.format(i)
        path_adj_out = dir_path + 'Adjout_{}.txt'.format(i)
        path_parameters = dir_path + 'cparameters_{}.txt'.format(i)
        # os.system(dir_path + '/slurm.serjob ' + dir_path + './cwesis.exe {} {}'.format(name_Adj_in,name_Adj_out))
        os.system('{} {} {} {}'.format(program_path,path_adj_in,path_adj_out,path_parameters))


def job_to_cluster(foldername,parameters,Istar):
    # This function submit jobs to the cluster with the following program keys:
    # bd: creates a bimodal directed networks and find its mean time to extinction
    dir_path = '/home/eladko/cpp_we_sis_extinction_project/we_sis_networks_extinction_cpp/'
    slurm_path = dir_path +'slurm.serjob'
    program_path = dir_path +'cwesis.exe'
    os.mkdir(foldername)
    os.chdir(foldername)
    data_path = os.path.dirname(os.path.realpath(__file__)) +'/' + foldername +'/'
    N,sims,it,k,x,lam,jump,Num_inf,Alpha,number_of_networks,tau,eps_din,eps_dout,new_trajcetory_bin,prog,Beta_avg = parameters
    if prog == 'bd':
        # G = nx.complete_graph(N)
        d1_in, d1_out, d2_in, d2_out = int(int(k) * (1 - float(eps_din))), int(int(k) * (1 - float(eps_dout))), int(int(k) * (1 + float(eps_din))), int(
            int(k) * (1 + float(eps_dout)))
        Beta = float(Beta_avg) / (1 + float(eps_din) * float(eps_dout))  # This is so networks with different std will have the reproduction number
        parameters = np.array([N,sims,it,k,x,lam,jump,Alpha,Beta,number_of_networks,tau,Istar,new_trajcetory_bin,dir_path,prog,dir_path])
        np.save('parameters.npy',parameters)
    # dir_path = os.path.dirname(os.path.realpath(__file__))
    for i in range(int(number_of_networks)):
        if prog=='bd':
            G = rand_networks.random_bimodal_directed_graph(int(d1_in), int(d1_out), int(d2_in), int(d2_out), int(N))
            parameters = np.array([N,sims,it,k,x,lam,jump,Alpha,Beta,i,tau,Istar,new_trajcetory_bin,dir_path,prog,dir_path])
        else:
            G = rand_networks.configuration_model_directed_graph(prog, float(eps_din), float(eps_dout), int(k), int(N))
            k_avg_graph = np.mean([G.in_degree(n) for n in G.nodes()])
            Beta_graph = float(lam)/k_avg_graph
            eps_in_graph = np.std([G.in_degree(n) for n in G.nodes()])/k_avg_graph
            eps_out_graph = np.std([G.out_degree(n) for n in G.nodes()])/k_avg_graph
            Beta = Beta_graph / (1 + np.sign(float(eps_din))*eps_in_graph * np.sign(float(eps_dout))* eps_out_graph)
            parameters = np.array([N,sims,it,k,x,lam,jump,Alpha,Beta,i,tau,Istar,new_trajcetory_bin,dir_path,prog,dir_path])
            np.save('parameters_{}.npy'.format(i), parameters)
        infile = 'GNull_{}.pickle'.format(i)
        with open(infile,'wb') as f:
            pickle.dump(G,f,pickle.HIGHEST_PROTOCOL)
        # nx.write_gpickle(G, infile)
        export_network_to_csv(G, i)
        export_parameters_to_csv(parameters,i)
        path_adj_in = dir_path + 'Adjin_{}.txt'.format(i)
        path_adj_out = dir_path + 'Adjout_{}.txt'.format(i)
        path_parameters = dir_path + 'cparameters_{}.txt'.format(i)
        parameters_path ='{} {} {}'.format(path_adj_in,path_adj_out,path_parameters)
        os.system('{} {} {}'.format(slurm_path,program_path,parameters_path))
        # os.system('{} {} {} {}'.format(program_path,path_adj_in,path_adj_out,path_parameters))



if __name__ == '__main__':
    # Parameters for the network to work
    N = 100 # number of nodes
    lam = 1.6 # The reproduction number
    number_of_networks = 10
    sims = 200 # Number of simulations at each step
    # k = N # Average number of neighbors for each node
    k = 50 # Average number of neighbors for each node
    x = 0.2 # intial infection percentage
    Num_inf = int(x*N) # Number of initially infected nodes
    it = 70
    jump = 5
    Alpha = 1.0 # Recovery rate
    Beta_avg = Alpha * lam / k # Infection rate for each node
    eps_din,eps_dout = 0.1,0.1 # The normalized std (second moment divided by the first) of the network
    # G = nx.random_regular_graph(k,N) # Creates a random graphs with k number of neighbors
    relaxation_time  = 3
    # tau = 1/(Num_inf*Alpha+N*Beta*k)
    tau = 1.0
    new_trajcetory_bin = 5
    prog = 'bd'
    dir_path = os.path.dirname(os.path.realpath(__file__))
    dir_path = '/home/elad/we_sis_networks_extinction_cpp/python'
    parameters = np.array([N,sims,it,k,x,lam,jump,Num_inf,Alpha,number_of_networks,tau,eps_din,eps_dout,new_trajcetory_bin,prog,Beta_avg])
    graphname  = 'GNull'
    foldername = 'prog_{}_k_{}_R_{}_tau_{}_it_{}_jump_{}_new_trajcetory_bin_{}_sims_{}_net_{}_epsin_{}_epsout_{}'.format(prog,k,lam,tau,it,jump,new_trajcetory_bin,sims,number_of_networks,eps_din,eps_dout)
    y1star=(-2*eps_din*(1 + eps_dout*eps_din)+ lam*(-1 + eps_din)*(1 + (-1 + 2*eps_dout)*eps_din)+ np.sqrt(lam**2 +eps_din*(4*eps_din +lam**2*eps_din*(-2 +eps_din**2) +4*eps_dout*(lam -(-2 + lam)*eps_din**2) +4*eps_dout**2*eps_din*(lam -(-1 + lam)*eps_din**2))))/(4*lam*(-1 +eps_dout)*(-1 +eps_din)*eps_din)
    y2star=(lam + eps_din*(-2 + 2*lam +lam*eps_din+ 2*eps_dout*(lam +(-1 + lam)*eps_din)) -np.sqrt(lam**2 +eps_din*(4*eps_din +lam**2*eps_din*(-2 +eps_din**2) +4*eps_dout*(lam -(-2 + lam)*eps_din**2) +4*eps_dout**2*eps_din*(lam -(-1 + lam)*eps_din**2))))/(4*lam*(1 +eps_dout)*eps_din*(1 + eps_din))
    Istar = (y1star +y2star)*N
    # Istar = (1 - 1/lam) * N


    # What's the job to run either on the cluster or on the laptop
    job_to_cluster(foldername,parameters,Istar)
    # act_as_main(foldername,parameters,Istar)

