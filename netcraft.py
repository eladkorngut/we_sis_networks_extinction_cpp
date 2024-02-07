import numpy as np
import rand_networks
import sys
import csv
import pickle
import networkx as nx


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


def create_network(file_name):
    parameters = np.load(file_name)
    N,sims,it,k,x,lam,jump,Num_inf,Alpha,Beta_avg,net_num,tau,new_trajcetory_bin,network_type,eps_din,eps_dout,Istar = parameters
    if network_type == 'bd':
        # G = nx.complete_graph(N)
        d1_in, d1_out, d2_in, d2_out = int(int(k) * (1 - float(eps_din))), int(int(k) * (1 - float(eps_dout))),\
                                       int(int(k) * (1 + float(eps_din))), int(int(k) * (1 + float(eps_dout)))
        Beta = float(Beta_avg) / (1 + float(eps_din) * float(eps_dout))  # This is so networks with different std will have the reproduction number
        parameters = np.array([N,sims,it,k,x,lam,jump,Alpha,Beta,net_num,tau,Istar,new_trajcetory_bin,network_type,eps_din,eps_dout])
        G = rand_networks.random_bimodal_directed_graph(int(d1_in), int(d1_out), int(d2_in), int(d2_out), int(N))
    elif network_type=='h':
        G = nx.random_regular_graph(int(k), int(N))
        parameters = np.array([N, sims, it, k, x, lam, jump, Alpha, Beta_avg, net_num, tau, Istar, new_trajcetory_bin, network_type, eps_din,eps_dout])  # Creates a random graphs with k number of neighbors
    else:
        G = rand_networks.configuration_model_directed_graph(network_type, float(eps_din), float(eps_dout), int(k),int(N))
        k_avg_graph = np.mean([G.in_degree(n) for n in G.nodes()])
        Beta_graph = float(lam) / k_avg_graph
        eps_in_graph = np.std([G.in_degree(n) for n in G.nodes()]) / k_avg_graph
        eps_out_graph = np.std([G.out_degree(n) for n in G.nodes()]) / k_avg_graph
        Beta = Beta_graph / (1 + np.sign(float(eps_din)) * eps_in_graph * np.sign(float(eps_dout)) * eps_out_graph)
        parameters = np.array([N, sims, it, k, x, lam, jump, Alpha, Beta, net_num, tau, Istar, new_trajcetory_bin,
                               network_type, eps_in_graph,eps_out_graph])
    np.save(file_name, parameters)
    infile = 'GNull_{}.pickle'.format(net_num)
    with open(infile, 'wb') as f:
        pickle.dump(G, f, pickle.HIGHEST_PROTOCOL)
    export_network_to_csv(G, net_num)
    export_parameters_to_csv(parameters,net_num)


if __name__ == '__main__':
    file_name = sys.argv[1]
    create_network(file_name)
