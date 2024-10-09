import numpy as np
import os
import rand_networks
import csv
import pickle
import networkx as nx
from scipy.stats import skew
import argparse
from scipy.sparse.linalg import eigsh

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

    # Check if the graph is directed
    is_directed = G.is_directed()

    with open('Adjin_{}.txt'.format(netname), 'w', newline='') as incoming_file:
        # Create a CSV writer
        incoming_writer = csv.writer(incoming_file)
        # Iterate over all nodes in the graph
        for node in np.sort(G):
            # Get the incoming neighbors of the current node
            if is_directed:
                incoming_neighbors = list(G.predecessors(node))
                # Get the degree of the current node
                degree = G.in_degree[node]
            else:
                incoming_neighbors = list(G.neighbors(node))  # All neighbors for undirected graph
                degree = G.degree[node]
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
            if is_directed:
                outgoing_neighbors = list(G.successors(node))
                # Get the degree of the current node
                degree = G.out_degree[node]
            else:
                outgoing_neighbors = list(G.neighbors(node))  # All neighbors for undirected graph
                degree = G.degree[node]
            # Write a row to the CSV file for the current node
            joint = np.concatenate(([degree],outgoing_neighbors),axis=0)
            outgoing_writer.writerow(joint)


def job_to_cluster(foldername,parameters,Istar,error_graphs,a):
    # This function submit jobs to the cluster with the following program keys:
    # bd: creates a bimodal skewed networks and find its mean time to extinction
    G, graph_degrees= 0,0
    dir_path = os.path.dirname(os.path.realpath(__file__))
    passed_error = False
    slurm_path = dir_path +'/slurm.serjob'
    program_path = dir_path +'/cwesis.exe'
    os.mkdir(foldername)
    os.chdir(foldername)
    data_path = os.getcwd() +'/'
    N, sims, it, k, x, lam, jump, Num_inf, Alpha, number_of_networks, tau, eps_din, eps_dout, new_trajcetory_bin, prog, Beta_avg, error_graphs_parameter,skewness = parameters
    N, sims, it, k, x, lam, jump, Num_inf, Alpha, number_of_networks, tau, eps_din, eps_dout, new_trajcetory_bin, prog, Beta_avg,skewness = \
        int(N), int(sims), int(it), float(k), float(x), float(lam), float(jump), int(Num_inf), float(Alpha), int(number_of_networks),\
        float(tau), float(eps_din), float(eps_dout),int(new_trajcetory_bin), prog, float(Beta_avg), float(skewness)
    for i in range(int(number_of_networks)):
        if not error_graphs and not passed_error:
            G,graph_degrees = rand_networks.configuration_model_undirected_graph_mulit_type(float(k),float(eps_din),int(N),prog,skewness)
            passed_error = True if error_graphs else False
            k_avg_graph,graph_std,graph_skewness,graph_median = np.mean(graph_degrees),np.std(graph_degrees),skew(graph_degrees),np.median(graph_degrees)
            second_moment,third_moment = np.mean((graph_degrees)**2),np.mean((graph_degrees)**3)
            eps_graph = graph_std / k_avg_graph
            N = len(G)
            largest_eigenvalue,largest_eigen_vector = eigsh(nx.adjacency_matrix(G).astype(float), k=1, which='LA', return_eigenvectors=True)
            Beta = float(lam) / largest_eigenvalue[0]
            graph_correlation = nx.degree_assortativity_coefficient(G)
            parameters = np.array(
                [N, sims, it, k_avg_graph, x, lam, jump, Alpha, Beta, i, tau, Istar, new_trajcetory_bin, dir_path,
                 prog, eps_graph, eps_graph, graph_std, graph_skewness, third_moment, second_moment,graph_correlation,graph_median])
        np.save('parameters_{}.npy'.format(i), parameters)
        np.save('largest_eigen_vector_{}.npy'.format(i), largest_eigen_vector)
        np.save('largest_eigenvalue_{}.npy'.format(i), largest_eigenvalue[0])
        infile = 'GNull_{}.pickle'.format(i)
        with open(infile,'wb') as f:
            pickle.dump(G,f,pickle.HIGHEST_PROTOCOL)
        export_network_to_csv(G, i)
        export_parameters_to_csv(parameters,i)
        path_adj_in = data_path + 'Adjin_{}.txt'.format(i)
        path_adj_out = data_path + 'Adjout_{}.txt'.format(i)
        path_parameters = data_path + 'cparameters_{}.txt'.format(i)
        parameters_path ='{} {} {}'.format(path_adj_in,path_adj_out,path_parameters)
        os.system('{} {} {}'.format(slurm_path,program_path,parameters_path))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process network and WE method parameters.")

    # Parameters for the network
    parser.add_argument('--N', type=int, help='Number of nodes')
    parser.add_argument('--prog', type=str, help='Program')
    parser.add_argument('--lam', type=float, help='The reproduction number')
    parser.add_argument('--eps_din', type=float, help='The normalized std (second moment divided by the first) of the in-degree distribution')
    parser.add_argument('--eps_dout', type=float, help='The normalized std (second moment divided by the first) of the out-degree distribution')
    parser.add_argument('--skewness', type=float, help='skewness parameter')
    parser.add_argument('--number_of_networks', type=int, help='Number of networks')
    parser.add_argument('--k', type=int, help='Average number of neighbors for each node')
    parser.add_argument('--error_graphs', action='store_true', help='Flag for error graphs')
    parser.add_argument('--a', type=float, help='Power-law distrbution first moment')

    # Parameters for the WE method
    parser.add_argument('--sims', type=int, help='Number of simulations at each bin')
    parser.add_argument('--tau', type=float, help='Tau parameter')
    parser.add_argument('--it', type=int, help='Number of iterations')
    parser.add_argument('--jump', type=int, help='Jump parameter')
    parser.add_argument('--new_trajectory_bin', type=int, help='New trajectory bin')

    # Parameters that don't get changed
    parser.add_argument('--relaxation_time', type=int, help='Relaxation time')
    parser.add_argument('--x', type=float, help='Initial infection percentage')
    parser.add_argument('--Alpha', type=float, help='Recovery rate')

    args = parser.parse_args()

    # Default parameters
    N = 6000 if args.N is None else args.N
    prog = 'pl' if args.prog is None else args.prog
    lam = 1.3 if args.lam is None else args.lam
    eps_din = 0.5 if args.eps_din is None else args.eps_din
    eps_dout = 0.5 if args.eps_dout is None else args.eps_dout
    skewness = 0.3 if args.skewness is None else args.skewness
    number_of_networks = 1 if args.number_of_networks is None else args.number_of_networks
    k = 50 if args.k is None else args.k
    a = 0.5 if args.a is None else args.a
    error_graphs = args.error_graphs

    sims = 500 if args.sims is None else args.sims
    tau = 1.0 if args.tau is None else args.tau
    it = 70 if args.it is None else args.it
    jump = 1 if args.jump is None else args.jump
    new_trajectory_bin = 2 if args.new_trajectory_bin is None else args.new_trajectory_bin

    relaxation_time = 20 if args.relaxation_time is None else args.relaxation_time
    x = 0.2 if args.x is None else args.x
    Num_inf = int(x * N)
    Alpha = 1.0 if args.Alpha is None else args.Alpha
    Beta_avg = Alpha * lam / k
    # run_mc_simulationtion = True


    parameters = np.array([N, sims, it, k, x, lam, jump, Num_inf, Alpha, number_of_networks, tau, eps_din, eps_dout, new_trajectory_bin, prog, Beta_avg, error_graphs, skewness])
    graphname = 'GNull'
    foldername = f'prog_{prog}_N{N}_k_{k}_R_{lam}_tau_{tau}_it_{it}_jump_{jump}_new_trajectory_bin_{new_trajectory_bin}_' \
                 f'sims_{sims}_net_{number_of_networks}_epsin_{eps_din}_epsout_{eps_dout}_skewness_{skewness}_a_{a}_err_{error_graphs}'
    Istar = (1 - 1/lam) * N

    job_to_cluster(foldername, parameters, Istar, error_graphs,a)