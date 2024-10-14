import csv

import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np
from collections import defaultdict
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from itertools import combinations
from itertools import zip_longest
from itertools import chain
from scipy.stats import skew


import numpy.random
from scipy.stats import norm
from scipy.stats import gamma
from scipy.stats import uniform
from scipy.stats import rv_discrete
import json


def draw_basic_nx_g(G):
    fig = plt.gcf()
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos,
                           node_color='r',
                           node_size=500,
                           alpha=1)
    nx.draw_networkx_edges(G, pos, width=1.0, alpha=1.0,arrowsize=20)
    plt.tight_layout()
    labels = nx.draw_networkx_labels(G, pos=pos,font_color='w')
    fig.savefig('network_strcture.png', dpi=100)
    plt.show()
    return pos


def draw_infections_nx_g(G,pos=None,headline=None,highest_rate=None):
    fig = plt.gcf()
    infected_list = []
    disinfected_list=[]
    for key, value in nx.get_node_attributes(G, 'infected').items():
        if value == True:
            infected_list.append(key)
        else:
            disinfected_list.append(key)
    if pos==None:
        pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos,
                           nodelist=infected_list,
                           node_color='r',
                           node_size=500,
                           alpha=1)
    nx.draw_networkx_nodes(G, pos,
                           nodelist=disinfected_list,
                           node_color='b',
                           node_size=500,
                           alpha=1)
    if highest_rate is not None:
        nx.draw_networkx_nodes(G, pos,
                               nodelist=highest_rate,
                               node_color='g',
                               node_size=500,
                               alpha=1)
    nx.draw_networkx_edges(G, pos, width=1.0, alpha=1.0)
    plt.tight_layout()
    labels = nx.draw_networkx_labels(G, pos=pos,font_color='w')
    fig.savefig('network_infection'+' '+str(headline)+'.png', dpi=100)
    plt.show()


# def draw_rate_nx_g(G):
#     Rate=[]
#     for n in range(G.number_of_nodes()):
#         Rate.append(G.nodes[n]['rate'])
#     Rate=np.dot(Rate,1/np.max(Rate))
#     Rate=np.array(Rate)
#     fig = plt.gcf()
#     pos = nx.spring_layout(G)
#     nx.draw_networkx_nodes(G, pos,
#                            cmap=plt.cm.Blues,
#                            node_size=500,
#                            alpha=1)
#     nx.draw_networkx_edges(G, pos, width=1.0, alpha=1.0)
#     plt.tight_layout()
#     labels = nx.draw_networkx_labels(G, pos=pos,font_color='w')
#     fig.savefig('/home/elad/Michael Assaf/binomal/photos/network_strcture.png', dpi=100)
#     plt.show()
#     return pos

def export_network(G):
    f = open('network_data.txt', "w+")
    for n in range(G.order()):
        f.write(str(G.degree(n)) + ' ')
        for neighbor in G[n]:
            f.write(str(neighbor) + ' ')
        f.write('\n')
    f.close()


def export_adj_network_cpp(G,infile):
    f = open(infile, "w+")
    for n in range(G.order()):
        f.write(str(G.degree(n)) + ',')
        for neighbor in G[n]:
            f.write(str(neighbor) + ',')
        f.write('\n')
    f.close()


def random_bimodal_graph(d1,d2, n, seed=None):
    """Return a random bimodal graph of n nodes each with degree d1 and d2.
        the degree distubition is half for each one of the two degrees

    The resulting graph G has no self-loops or parallel edges.

    Parameters
    ----------
    d1 : int
      Degre
    d2 : int
      Degree
    n : integer
      Number of nodes. The value of n*d must be even.
    seed : hashable object
        The seed for random number generator.

    Notes
    -----
    The nodes are numbered form 0 to n-1.

    Kim and Vu's paper [2]_ shows that this algorithm samples in an
    asymptotically uniform way from the space of random graphs when
    d = O(n**(1/3-epsilon)).

    References
    ----------
    .. [1] A. Steger and N. Wormald,
       Generating random regular graphs quickly,
       Probability and Computing 8 (1999), 377-396, 1999.
       http://citeseer.ist.psu.edu/steger99generating.html

    .. [2] Jeong Han Kim and Van H. Vu,
       Generating random regular graphs,
       Proceedings of the thirty-fifth ACM symposium on Theory of computing,
       San Diego, CA, USA, pp 213--222, 2003.
       http://portal.acm.org/citation.cfm?id=780542.780576

       The multiplcation of n*d1 and n*d2 needs to be even. For the safe side choose even number for both nodes and degree
    """
    if (n * d1) % 2 != 0 or (n * d2) % 2!=0:
        raise nx.NetworkXError("n * d must be even")

    if not (0 <= d1 < n or 0 <= d2 < n):
        raise nx.NetworkXError("the 0 <= d < n inequality must be satisfied")

    if d1 == 0 or d2 == 0:
        return nx.empty_graph(n)

    if seed is not None:
        random.seed(seed)

    def _suitable(edges, potential_edges):
    # Helper subroutine to check if there are suitable edges remaining
    # If False, the generation of the graph has failed
        if not potential_edges:
            return True
        for s1 in potential_edges:
            for s2 in potential_edges:
                # Two iterators on the same dictionary are guaranteed
                # to visit it in the same order if there are no
                # intervening modifications.
                if s1 == s2:
                    # Only need to consider s1-s2 pair one time
                    break
                if s1 > s2:
                    s1, s2 = s2, s1
                if (s1, s2) not in edges:
                    return True
        return False

    def _try_creation(n):
        # Attempt to create an edge set

        edges = set()
        stub_d1=random.sample(list(range(n)),int(n/2))
        stub_d2= [item for item in list(range(n)) if item not in stub_d1]
        stub_d1=stub_d1 * d1
        stub_d2=stub_d2 * d2
        stubs = stub_d1+stub_d2

        while stubs:
            potential_edges = defaultdict(lambda: 0)
            random.shuffle(stubs)
            stubiter = iter(stubs)
            for s1, s2 in zip(stubiter, stubiter):
                if s1 > s2:
                    s1, s2 = s2, s1
                if s1 != s2 and ((s1, s2) not in edges):
                    edges.add((s1, s2))
                else:
                    potential_edges[s1] += 1
                    potential_edges[s2] += 1

            if not _suitable(edges, potential_edges):
                return None # failed to find suitable edge set

            stubs = [node for node, potential in potential_edges.items()
                     for _ in range(potential)]
        return edges

    # Even though a suitable edge set exists,
    # the generation of such a set is not guaranteed.
    # Try repeatedly to find one.
    edges = _try_creation(n)
    while edges is None:
        edges = _try_creation(n)

    G = nx.Graph()
    G.name = "random_bimodal_graph(%s, %s)" % (d1, n)
    G.add_edges_from(edges)

    return G


def random_bimodal_graph_skewed(d1,d2, n, fraction,seed=None):
    """Return a random bimodal graph of n nodes each with degree d1 and d2.
        the degree distubition is half for each one of the two degrees

    The resulting graph G has no self-loops or parallel edges.

    Parameters
    ----------
    d1 : int
      Degre
    d2 : int
      Degree
    n : integer
      Number of nodes. The value of n*d must be even.
    seed : hashable object
        The seed for random number generator.

    Notes
    -----
    The nodes are numbered form 0 to n-1.

    Kim and Vu's paper [2]_ shows that this algorithm samples in an
    asymptotically uniform way from the space of random graphs when
    d = O(n**(1/3-epsilon)).

    References
    ----------
    .. [1] A. Steger and N. Wormald,
       Generating random regular graphs quickly,
       Probability and Computing 8 (1999), 377-396, 1999.
       http://citeseer.ist.psu.edu/steger99generating.html

    .. [2] Jeong Han Kim and Van H. Vu,
       Generating random regular graphs,
       Proceedings of the thirty-fifth ACM symposium on Theory of computing,
       San Diego, CA, USA, pp 213--222, 2003.
       http://portal.acm.org/citation.cfm?id=780542.780576

       The multiplcation of n*d1 and n*d2 needs to be even. For the safe side choose even number for both nodes and degree
    """
    if (n * d1) % 2 != 0 or (n * d2) % 2!=0:
        raise nx.NetworkXError("n * d must be even")

    if not (0 <= d1 < n or 0 <= d2 < n):
        raise nx.NetworkXError("the 0 <= d < n inequality must be satisfied")

    if d1 == 0 or d2 == 0:
        return nx.empty_graph(n)

    if seed is not None:
        random.seed(seed)

    def _suitable(edges, potential_edges):
    # Helper subroutine to check if there are suitable edges remaining
    # If False, the generation of the graph has failed
        if not potential_edges:
            return True
        for s1 in potential_edges:
            for s2 in potential_edges:
                # Two iterators on the same dictionary are guaranteed
                # to visit it in the same order if there are no
                # intervening modifications.
                if s1 == s2:
                    # Only need to consider s1-s2 pair one time
                    break
                if s1 > s2:
                    s1, s2 = s2, s1
                if (s1, s2) not in edges:
                    return True
        return False

    def _try_creation(n):
        # Attempt to create an edge set

        edges = set()
        stub_d1 = random.sample(list(range(n)),int(n*fraction))
        stub_d2= [item for item in list(range(n)) if item not in stub_d1]
        stub_d1=stub_d1 * d1
        stub_d2=stub_d2 * d2
        stubs = stub_d1+stub_d2

        while stubs:
            potential_edges = defaultdict(lambda: 0)
            random.shuffle(stubs)
            stubiter = iter(stubs)
            for s1, s2 in zip(stubiter, stubiter):
                if s1 > s2:
                    s1, s2 = s2, s1
                if s1 != s2 and ((s1, s2) not in edges):
                    edges.add((s1, s2))
                else:
                    potential_edges[s1] += 1
                    potential_edges[s2] += 1

            if not _suitable(edges, potential_edges):
                return None # failed to find suitable edge set

            stubs = [node for node, potential in potential_edges.items()
                     for _ in range(potential)]
        return edges

    # Even though a suitable edge set exists,
    # the generation of such a set is not guaranteed.
    # Try repeatedly to find one.
    edges = _try_creation(n)
    while edges is None:
        edges = _try_creation(n)

    G = nx.Graph()
    G.name = "random_bimodal_graph(%s, %s)" % (d1, n)
    G.add_edges_from(edges)

    return G



def random_bimodal_directed_graph(d1_in,d1_out,d2_in,d2_out, n, seed=None):
    """Return a random bimodal directed graph of n nodes each with degree d1 and d2.
        the degree distubition is half for each one of the two degrees

    The resulting graph G has no self-loops or parallel edges.

    Parameters
    ----------
    d1 : int
      Degre
    d2 : int
      Degree
    n : integer
      Number of nodes. The value of n*d must be even.
    seed : hashable object
        The seed for random number generator.

    Notes
    -----
    The nodes are numbered form 0 to n-1.

    Kim and Vu's paper [2]_ shows that this algorithm samples in an
    asymptotically uniform way from the space of random graphs when
    d = O(n**(1/3-epsilon)).

    References
    ----------
    .. [1] A. Steger and N. Wormald,
       Generating random regular graphs quickly,
       Probability and Computing 8 (1999), 377-396, 1999.
       http://citeseer.ist.psu.edu/steger99generating.html

    .. [2] Jeong Han Kim and Van H. Vu,
       Generating random regular graphs,
       Proceedings of the thirty-fifth ACM symposium on Theory of computing,
       San Diego, CA, USA, pp 213--222, 2003.
       http://portal.acm.org/citation.cfm?id=780542.780576

       The multiplcation of n*d1 and n*d2 needs to be even. For the safe side choose even number for both nodes and degree
    """
    if (n * d1_in) % 2 != 0 or (n * d2_in) % 2!=0:
        raise nx.NetworkXError("n * d must be even")

    if not (0 <= d1_in < n or 0 <= d2_in < n):
        raise nx.NetworkXError("the 0 <= d < n inequality must be satisfied")

    if not (d1_in+d2_in==d1_out+d2_out):
        raise nx.NetworkXError("the d1_in+d2_in==d1_out+d2_out equality must be satisfied")

    if seed is not None:
        random.seed(seed)

    def _suitable(edges, potential_edges_in,potential_edges_out):
    # Helper subroutine to check if there are suitable edges remaining
    # If False, the generation of the graph has failed
        if not (potential_edges_in and potential_edges_out):
            return True
        for s1 in potential_edges_in:
            for s2 in potential_edges_out:
                # Two iterators on the same dictionary are guaranteed
                # to visit it in the same order if there are no
                # intervening modifications.
                if s1 == s2:
                    # Only need to consider s1-s2 pair one time
                    break
                if (s1, s2) not in edges:
                    return True
        return False

    def _try_creation(n):
        # Attempt to create an edge set

        edges = set()
        stub_d1_in = random.sample(list(range(n)),int(n/2))
        stub_d2_in = [item for item in list(range(n)) if item not in stub_d1_in]
        stub_d1_in =stub_d1_in * d1_in
        stub_d2_in = stub_d2_in * d2_in
        stubs_in = stub_d1_in+stub_d2_in

        stub_d1_out = [item for item in list(range(n)) if item in stub_d1_in]
        stub_d2_out = [item for item in list(range(n)) if item in stub_d2_in]
        stub_d1_out = stub_d1_out * d1_out
        stub_d2_out = stub_d2_out * d2_out
        stubs_out = stub_d1_out + stub_d2_out


        while stubs_in and stubs_out:
            potential_edges_in = defaultdict(lambda: 0)
            potential_edges_out = defaultdict(lambda: 0)
            random.shuffle(stubs_in)
            stubiter_in = iter(stubs_in)
            random.shuffle(stubs_out)
            stubiter_out = iter(stubs_out)
            for s1, s2 in zip(stubiter_in, stubiter_out):
                if s1 != s2 and ((s1, s2) not in edges):
                    edges.add((s1, s2))
                else:
                    potential_edges_in[s1] += 1
                    potential_edges_out[s2] += 1

            if not _suitable(edges, potential_edges_in,potential_edges_out):
                return None # failed to find suitable edge set

            stubs_in = [node for node, potential in potential_edges_in.items()
                     for _ in range(potential)]
            stubs_out = [node for node, potential in potential_edges_out.items()
                     for _ in range(potential)]
        return edges

    # Even though a suitable edge set exists,
    # the generation of such a set is not guaranteed.
    # Try repeatedly to find one.
    edges = _try_creation(n)
    while edges is None:
        edges = _try_creation(n)

    G = nx.DiGraph()
    G.name = "random_bimodal_graph(%s, %s)" % (d1_in, n)
    G.add_edges_from(edges)
    return G

def create_exact_gauss_dis(epsilon,avg_degree,N):
    low, high,normalization = (avg_degree * (1 - 4 * epsilon)).astype(int), (avg_degree * (1 + 4 * epsilon)).astype(int),N
    possible_degrees = np.arange(low, high)
    deg_dist = (norm.pdf(possible_degrees, loc=avg_degree, scale=avg_degree * epsilon) * normalization).astype(int)
    deg_for_algo = list(chaini([n] * d for n, d in zip(possible_degrees, deg_dist)))
    while np.size(deg_for_algo)<N:
        normalization=normalization+1
        deg_dist = (norm.pdf(np.arange(low, high), loc=avg_degree, scale=avg_degree * epsilon) * normalization).astype(int)
        deg_for_algo = list(chaini([n] * d for n, d in zip(possible_degrees, deg_dist)))
    deg_dist, possible_degrees = deg_dist[deg_dist != 0], possible_degrees[deg_dist != 0]
    for i in range(np.size(deg_for_algo)-N): deg_dist[i] =deg_dist[i] - 1
    # deg_dist[0] = deg_dist[0] - 1 if np.size(deg_for_algo) != N else deg_dist[0]
    deg_for_algo = list(chaini([n] * d for n, d in zip(possible_degrees, deg_dist)))
    return np.array(deg_for_algo)


def create_deg_distubtion_conf_model(epsilon,avg_degree,N,dist_type):

    def create_stubs_degree(epsilon,avg_degree,dist_type,normalization):
        if dist_type=='gauss_c':
            low, high = (avg_degree * (1 - 4 * epsilon)).astype(int), (
                        avg_degree * (1 + 4 * epsilon)).astype(int)
            possible_degrees = np.arange(low, high)
            deg_dist = (norm.pdf(possible_degrees, loc=avg_degree, scale=avg_degree * epsilon) * normalization).astype(int)
        elif dist_type == 'bimodal_c':
            possible_degrees = np.array([(avg_degree * (1 - epsilon)).astype(int), (avg_degree * (1 + epsilon)).astype(int)])
            deg_dist = np.array([N/2,N/2]).astype(int)
        elif dist_type == 'gamma_c':
            low, high = (avg_degree * (1 - 100 * epsilon)).astype(int), (
                        avg_degree * (1 + 100 * epsilon)).astype(int)
            possible_degrees = np.arange(low, high)
            deg_dist = (gamma.pdf(possible_degrees, a=1/epsilon**2,scale=epsilon**2*avg_degree)*normalization).astype(int)
        elif dist_type == 'uniform_c':
            low,high = (avg_degree * (1 - 100 * epsilon)).astype(int), (
                        avg_degree * (1 + 100 * epsilon)).astype(int)
            possible_degrees = np.arange(low, high)
            deg_dist = (uniform.pdf(possible_degrees, loc=int(avg_degree*(1-epsilon)),scale=2*np.sqrt(3)*epsilon*avg_degree)*normalization).astype(int)
        return list(chaini([n] * d for n, d in zip(possible_degrees, deg_dist))),deg_dist,possible_degrees

    def trim_edges_to_fit_N_nodes(deg_for_algo,deg_dist,possible_degrees,N):
        deg_dist, possible_degrees = deg_dist[deg_dist != 0], possible_degrees[deg_dist != 0]
        reduce_by_one =  np.random.uniform(0, np.size(deg_dist)-1, np.size(deg_for_algo) - N).astype(int)
        for i in reduce_by_one: deg_dist[i]=deg_dist[i]-1
        # for i in range(np.size(deg_for_algo) - N): deg_dist[i] = deg_dist[i] - 1
        return deg_dist,possible_degrees

    normalization=N
    deg_for_algo,deg_dist,possible_degrees=create_stubs_degree(epsilon,avg_degree,dist_type,normalization)
    while np.size(deg_for_algo)<=N:
        normalization=normalization+1
        deg_for_algo,deg_dist,possible_degrees = create_stubs_degree(epsilon,avg_degree,dist_type,normalization)
    deg_dist, possible_degrees = trim_edges_to_fit_N_nodes(deg_for_algo,deg_dist,possible_degrees,N)
    deg_for_algo = list(chaini([n] * d for n, d in zip(possible_degrees, deg_dist)))
    return np.array(deg_for_algo)




def create_random_degree_sequence(name,epsilon,avg_degree,N):
    if epsilon==0.0:
        return np.ones(N).astype(int)*int(avg_degree)
    if name=='gauss':
        return np.random.normal(avg_degree, epsilon * avg_degree, N)
    elif name=='uniform':
        return np.random.uniform(avg_degree*(1-np.sqrt(3)*epsilon),avg_degree*(1+np.sqrt(3)*epsilon),N).astype(int)
    elif name=='gamma':
        return np.random.gamma(1/epsilon**2,avg_degree*epsilon**2,N)
    elif name=='gauss_e':
        return create_exact_gauss_dis(epsilon,avg_degree,N)
    elif name=='gauss_c' or name=='bimodal_c' or name=='gamma_c' or name=='uniform_c':
        return create_deg_distubtion_conf_model(epsilon,avg_degree,N,name)


def equate_degree_sequence(din,dout,N,avg_degree,eps_in=0.0,eps_out=0.0):
    din[din < 0] = 0
    dout[dout < 0] = 0
    return_order_din_dout=True
    if eps_in==0.0:
        for i in range(N*avg_degree-np.sum(dout)):
            node = int(N * np.random.random())
            dout[node] = dout[node] + 1
        return din,dout
    if eps_out==0.0:
        for i in range(N*avg_degree-np.sum(din)):
            node = int(N * np.random.random())
            din[node] = din[node] + 1
        return din,dout
    if np.sum(din)>np.sum(dout):
        high, low = din, dout
    else:
        return_order_din_dout = False
        high, low = dout, din
    prob=np.random.uniform(0,1,np.sum(high)-np.sum(low))
    for p in prob:
        if p>0.5:
            node=int(N*np.random.random())
            while high[node]==0:
                node = int(N * np.random.random())
            high[node]=high[node]-1
        else:
            node=int(N*np.random.random())
            while low[node]>=N:
                node = int(N * np.random.random())
            low[node]=low[node]+1
    if return_order_din_dout: return high,low
    return low,high


chaini = chain.from_iterable

def _to_stublist_simple_graph(degree_sequence):
    """Returns a list of degree-repeated node numbers.

    ``degree_sequence`` is a list of nonnegative integers representing
    the degrees of nodes in a graph.

    This function returns a list of node numbers with multiplicities
    according to the given degree sequence. For example, if the first
    element of ``degree_sequence`` is ``3``, then the first node number,
    ``0``, will appear at the head of the returned list three times. The
    node numbers are assumed to be the numbers zero through
    ``len(degree_sequence) - 1``.

    Examples
    --------

    # >>> degree_sequence = [1, 2, 3]
    # >>> _to_stublist(degree_sequence)
    [0, 1, 1, 2, 2, 2]

    If a zero appears in the sequence, that means the node exists but
    has degree zero, so that number will be skipped in the returned
    list::

    # >>> degree_sequence = [2, 0, 1]
    # >>> _to_stublist(degree_sequence)
    [0, 0, 2]

    """
    return list(chaini([n] * d for n, d in enumerate(degree_sequence)))


def _configuration_model_simple_graph(
    create_using,name,epsilon_in,epsilon_out,avg_degree,n):
    """Helper function for generating either undirected or directed
    configuration model graphs.

    ``deg_sequence`` is a list of nonnegative integers representing the
    degree of the node whose label is the index of the list element.

    ``create_using`` see :func:`~networkx.empty_graph`.

    ``directed`` and ``in_deg_sequence`` are required if you want the
    returned graph to be generated using the directed configuration
    model algorithm. If ``directed`` is ``False``, then ``deg_sequence``
    is interpreted as the degree sequence of an undirected graph and
    ``in_deg_sequence`` is ignored. Otherwise, if ``directed`` is
    ``True``, then ``deg_sequence`` is interpreted as the out-degree
    sequence and ``in_deg_sequence`` as the in-degree sequence of a
    directed graph.

    .. note::

       ``deg_sequence`` and ``in_deg_sequence`` need not be the same
       length.

    ``seed`` is a random.Random or numpy.random.RandomState instance

    This function returns a graph, directed if and only if ``directed``
    is ``True``, generated according to the configuration model
    algorithm. For more information on the algorithm, see the
    :func:`configuration_model` or :func:`directed_configuration_model`
    functions.

    """
    G = nx.empty_graph(n, create_using)
    # If empty, return the null graph immediately.
    if n == 0:
        return G
    # Build a list of available degree-repeated nodes.  For example,
    # for degree sequence [3, 2, 1, 1, 1], the "stub list" is
    # initially [0, 0, 0, 1, 1, 2, 3, 4], that is, node 0 has degree
    # 3 and thus is repeated 3 times, etc.
    #
    # Also, shuffle the stub list in order to get a random sequence of
    # node pairs.
    def _suitalbe(edges,potential_edges_in,potential_edges_out):
        # Helper suroutine to check if there are suitable edges remaining
        # If false generation of graph failed
        if not (potential_edges_in and potential_edges_out):
            return True
        for s1 in potential_edges_in:
            for s2 in potential_edges_out:
                if s1==s2:
                    break
                if (s2,s1) not in edges:
                    return True
        return False
    def _try_creation(n,epsilon_in,epsilon_out,avg_degree):
        edges = set()
        din, dout = equate_degree_sequence(create_random_degree_sequence(name, np.abs(epsilon_in), avg_degree, n),
                                           create_random_degree_sequence(name, np.abs(epsilon_out), avg_degree, n), n,avg_degree,epsilon_in,epsilon_out)
        din = np.sort(din)
        dout = np.sort(dout) if epsilon_out * epsilon_in > 0 else np.sort(dout)[::-1]
        pairs = zip_longest(dout, din, fillvalue=0)
        # Unzip the list of pairs into a pair of lists.
        out_deg, in_deg = zip(*pairs)
        out_stublist = _to_stublist_simple_graph(out_deg)
        in_stublist = _to_stublist_simple_graph(in_deg)
        while in_stublist and out_stublist:
            potential_edge_in = defaultdict(lambda: 0)
            potential_edge_out = defaultdict(lambda: 0)
            np.random.shuffle(in_stublist)
            subiter_in = iter(in_stublist)
            np.random.shuffle(out_stublist)
            subiter_out = iter(out_stublist)
            for s1,s2 in zip(subiter_in,subiter_out):
                if s1!=s2 and ((s2,s1) not in edges):
                    edges.add((s2,s1))
                else:
                    potential_edge_in[s1] += 1
                    potential_edge_out[s2] += 1
            if not _suitalbe(edges,potential_edge_in,potential_edge_out):
                return None # failed to find suitable edge set
            in_stublist = [node for node, potential in potential_edge_in.items() for _ in range(potential)]
            out_stublist = [node for node,potential in potential_edge_out.items() for _ in range(potential)]
        return edges
    edges = _try_creation(n,epsilon_in,epsilon_out,avg_degree)
    while edges is None:
        edges = _try_creation(n,epsilon_in,epsilon_out,avg_degree)
    G =nx.DiGraph()
    G.name = "simple_directed_config_graph"
    G.add_edges_from(edges)
    return G



def configuration_model_directed_graph(
    name, epsilon_in,epsilon_out,avg_degree,n, create_using=None):
    """Returns a directed_random graph with the given degree sequences.

    The configuration model generates a random directed pseudograph
    by randomly assigning edges to match the given degree sequences.
    Self edges and loops are excluded by randomly choosing
    Parameters
    ----------
    in_degree_sequence :  list of nonnegative integers
       Each list entry corresponds to the in-degree of a node.
    out_degree_sequence :  list of nonnegative integers
       Each list entry corresponds to the out-degree of a node.
    create_using : NetworkX graph constructor, optional (default MultiDiGraph)
        Graph type to create. If graph instance, then cleared before populated.
    seed : integer, random_state, or None (default)
        Indicator of random number generation state.
        See :ref:`Randomness<randomness>`.

    Returns
    -------
    G : MultiDiGraph
        A graph with the specified degree sequences.
        Nodes are labeled starting at 0 with an index
        corresponding to the position in deg_sequence.

    Raises
    ------
    NetworkXError
        If the degree sequences do not have the same sum.

    See Also
    --------
    configuration_model

    Notes
    -----
    Algorithm as described by Newman [1]_.

    References
    ----------
    .. [1] Newman, M. E. J. and Strogatz, S. H. and Watts, D. J.
       Random graphs with arbitrary degree distributions and their applications
       Phys. Rev. E, 64, 026118 (2001)
    """

    if create_using is None:
        create_using = nx.MultiDiGraph

    G = _configuration_model_simple_graph(
        create_using,
        name,
        epsilon_in,
        epsilon_out,
        avg_degree,
        n
    )

    name = "directed configuration_model {} nodes {} edges"
    return G

def plot_configuration_model_powerlaw_multi(G_array,a_array,b,n,custom_dist,k):
    marker_type, line_style = ['o', 'v', '+'], ['-', '--', ':']
    for G, a, mk, ls in zip(G_array, a_array, marker_type, line_style):
        degree_sequence = sorted((d for n, d in G.degree()), reverse=True)
        hist, edges = np.histogram(degree_sequence, bins=np.arange(1, np.max(degree_sequence)))
        plt.plot(edges[:-1], hist, marker=mk, linestyle='none', label='a={}'.format(a))
        # x = np.arange(1, n)
        # pmf_values = n * custom_dist.pmf(x, a, b)

        # Filter values where pmf_values is greater than zero
        # valid_indices = pmf_values > 1
        # x_filtered = x[valid_indices]
        # pmf_values_filtered = pmf_values[valid_indices]

        # plt.plot(x_filtered, pmf_values_filtered, 'b-', label='PMF', linewidth=2)
    # Define the start and end points of the lines
    x_start, y_start = 1, 10000
    x_end, y_end = 700, 1

    x_start_2, y_start_2 = 3, 1100
    x_end_2, y_end_2 = 300, 1

    # Convert the start and end points to log scale
    log_x_start = np.log10(x_start)
    log_y_start = np.log10(y_start)
    log_x_end = np.log10(x_end)
    log_y_end = np.log10(y_end)

    log_x_start2 = np.log10(x_start_2)
    log_y_start2 = np.log10(y_start_2)
    log_x_end2 = np.log10(x_end_2)
    log_y_end2 = np.log10(y_end_2)

    # Fit a linear line to the log-transformed data for slope -2
    m2 = -1.5  # Slope of the log-transformed line for slope -2
    log_x2 = np.linspace(log_x_start2, log_x_end2, 100)
    log_y2 = m2 * (log_x2 - log_x_start2) + log_y_start2

    # Convert the log-scale data back to linear scale for slope -2
    x2 = 10 ** log_x2
    y2 = 10 ** log_y2

    # Plot the line for slope -2
    plt.plot(x2, y2, color='k', linestyle='-', label='Slope -1.5')

    # Fit a linear line to the log-transformed data for slope -1
    m1 = -1  # Slope of -1
    log_x1 = np.linspace(log_x_start, log_x_end, 100)
    log_y1 = m1 * (log_x1 - log_x_start) + log_y_start

    # Convert the log-scale data back to linear scale for slope -1
    x1 = 10 ** log_x1
    y1 = 10 ** log_y1

    # Plot the line for slope -1
    plt.plot(x1, y1, color='r', linestyle='--', label='Slope -1')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Degree (k)')
    plt.ylabel('P(k)')
    plt.title('Power-law Degree Distribution with N={}'.format(n))
    plt.legend()
    plt.grid(True)
    plt.savefig('power_law_graph_all_a.png', dpi=500)
    plt.show()

def plot_configuration_model_powerlaw(G,a,b,n,custom_dist,k):
    degree_sequence = sorted((d for n, d in G.degree()), reverse=True)
    hist, edges = np.histogram(degree_sequence, bins=np.arange(1, np.max(degree_sequence)))
    plt.plot(edges[:-1], hist, 'ro', label='Network')
    x = np.arange(1, n)
    pmf_values = n * custom_dist.pmf(x, a, b)

    # Filter values where pmf_values is greater than zero
    valid_indices = pmf_values > 1
    x_filtered = x[valid_indices]
    pmf_values_filtered = pmf_values[valid_indices]

    plt.plot(x_filtered, pmf_values_filtered, 'b-', label='PMF', linewidth=2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Degree (k)')
    plt.ylabel('P(k)')
    plt.title('Power-law Degree Distribution with k={}, a={}, b={} and N={}'.format(round(np.mean(degree_sequence), 2),
                                                                                    round(a, 2), round(b, 2), n))
    plt.legend()
    plt.grid(True)
    plt.savefig('power_law_graph_a{}.png'.format(round(a,1)), dpi=500)
    plt.show()

def plot_gamma_distribution(G,kavg,epsilon,n,net_type):
    if net_type=='bet':
        title='Beta'
    elif net_type=='gam':
        title='Gamma'
    elif net_type=='ig':
        title = 'Wald'
    elif net_type=='ln':
        title='Log-norm'
    degree_sequence = sorted((d for n, d in G.degree()), reverse=True)
    hist, edges = np.histogram(degree_sequence, bins=np.arange(1, np.max(degree_sequence)),density=False)
    plt.plot(edges[:-1], hist, 'ro', label='Network <k>={}'.format(np.mean(degree_sequence)))
    possible_degrees = range(0,max(degree_sequence))
    # gamma_theory = (gamma.pdf(possible_degrees, a=1 / epsilon ** 2, scale=epsilon ** 2 * kavg)*n)
    # plt.plot(possible_degrees,gamma_theory,'-b',label='Theory')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlabel('Degree (k)')
    plt.ylabel('P(k)')
    plt.title(f'N={len(degree_sequence)}, <k>={round(np.mean(degree_sequence),2)}, Median={np.median(degree_sequence)} skewness={round(skew(degree_sequence),2)} eps={round(np.std(degree_sequence)/np.mean(degree_sequence),2)}')
    # plt.title(r'Degree Distribution {} with k={}, $\epsilon$={} and N={}'.format(title,round(np.mean(degree_sequence),2), round(epsilon,2), n))
    plt.legend()
    plt.grid(True)
    plt.savefig('degree_dist_graph_{}_k{}_eps_{}_N{}.png'.format(net_type,round(kavg,1),round(epsilon,2),n), dpi=500)
    plt.show()

def find_a_binary_search(kavg, n, a):
    class CustomDistribution(rv_discrete):
        def _pmf(self, k, a, b):
            return b * a / (1 + b * k) ** (a + 1)
    b = 1/(kavg*(a-1))
    custom_dist =CustomDistribution()
    mid_mean = np.mean(custom_dist.rvs(a=a, b=b, size=n))
    high, low = a, a - 2
    while np.abs(mid_mean - kavg) / kavg > 0.01:
        mid = (low + high) / 2
        mid_mean = np.mean(custom_dist.rvs(a=mid, b=b, size=n))
        if mid_mean > kavg:
            low = mid
        else:
            high = mid
    return (low + high) / 2,b


def find_b_binary_search(kavg,n,a):
    class CustomDistribution(rv_discrete):
        def _pmf(self, k, a, b):
            return b * a / (1 + b * k) ** (a + 1)
    custom_dist =CustomDistribution()
    high, low = 10.0, 0.0
    mid = (low + high) / 2
    mid_mean = np.mean(custom_dist.rvs(a=a, b=mid, size=n))
    while np.abs(mid_mean - kavg) / kavg > 0.01:
        mid = (low + high) / 2
        mid_mean = np.mean(custom_dist.rvs(a=a, b=mid, size=n))
        if mid_mean > kavg:
            low = mid
        else:
            high = mid
    return a,(low + high) / 2

def configuration_model_powerlaw(a,b,n):
    class CustomDistribution(rv_discrete):
        def _pmf(self, k, a, b):
            return b * a / (1 + b * k) ** (a + 1)
    custom_dist =CustomDistribution()
    sequence = custom_dist.rvs(a=a,b=b,size=n)
    return sequence
    # if (np.sum(sequence)%2!=0):
    #     rand_bin = random.randint(1, n)
    #     sequence[rand_bin]=sequence[rand_bin]+1
    # G = nx.configuration_model(sequence)
    # G = nx.Graph(G)
    # G.remove_edges_from(nx.selfloop_edges(G))
    # return G


def configuration_model_powerlaw_degree(a,b,n):
    class CustomDistribution(rv_discrete):
        def _pmf(self, k, a, b):
            return b * a / (1 + b * k) ** (a + 1)
    custom_dist = CustomDistribution()
    sequence = custom_dist.rvs(a=a,b=b,size=n)
    if (np.sum(sequence)%2!=0):
        rand_bin = random.randint(1, n)
        sequence[rand_bin]=sequence[rand_bin]+1
    return G



def configuration_model_undirected_graph(epsilon,avg_degree,N):
    d=np.random.normal(avg_degree, epsilon * avg_degree, N).astype(int)
    if np.sum(d)%2!=0:
        d[int(N*np.random.random())]+=1
    G = nx.configuration_model(d.astype(int))
    G = nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    return G,a,b

def configuration_model_undirected_graph_exp(avg_degree,N):
    d=np.random.exponential(avg_degree, N).astype(int)
    if np.sum(d)%2!=0:
        d[int(N*np.random.random())]+=1
    G = nx.configuration_model(d.astype(int))
    G = nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    return G


def naive_config_gamma(kavg,epsilon,N):
    theta, shape, k_avg_graph = epsilon ** 2 * kavg, 1 / epsilon ** 2, 0.0
    d = numpy.random.default_rng().gamma(shape, theta, N).astype(int)
    if np.sum(d) % 2 != 0:
        d[int(len(d) * np.random.random())] += 1
    G = nx.configuration_model(d)
    G = nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    k_avg_graph = np.mean([G.degree(n) for n in G.nodes()])
    return G,k_avg_graph

def naive_config_multi(kavg,epsilon,N,net_type):
    if net_type == 'ig':
        wald_mu, wald_lambda = kavg, kavg / epsilon ** 2
        d = numpy.random.default_rng().wald(wald_mu, wald_lambda, N).astype(int)
    elif net_type == 'bet':
        alpha_beta_dist, beta_beta_dist = (N - kavg * (1 + epsilon ** 2)) / (N * epsilon ** 2), (
                    (kavg - N) * (kavg - N + kavg * epsilon ** 2)) / (kavg * N * epsilon ** 2)
        d = (numpy.random.default_rng().beta(alpha_beta_dist, beta_beta_dist, N) * N).astype(int)
    elif net_type == 'ln':
        mu_log_norm, sigma_log_norm = -(1 / 2) * np.log((1 + epsilon ** 2) / kavg ** 2), np.sqrt(
            2 * np.log(kavg) + np.log((1 + epsilon ** 2) / kavg ** 2))
        d = numpy.random.lognormal(mu_log_norm, sigma_log_norm, N).astype(int)
    elif net_type == 'gam':
        theta, shape, k_avg_graph = epsilon ** 2 * kavg, 1 / epsilon ** 2, 0.0
        d = numpy.random.default_rng().gamma(shape, theta, N).astype(int)
    if np.sum(d) % 2 != 0:
        d[int(len(d) * np.random.random())] += 1
    G = nx.configuration_model(d)
    G = nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    k_avg_graph = np.mean([G.degree(n) for n in G.nodes()])
    return G,k_avg_graph

def find_k_binary_search(kavg,epsilon,n):
    high, low = n, kavg
    mid = (low + high) / 2
    G_mean,mid_mean = naive_config_gamma(mid,epsilon,n)
    while np.abs(mid_mean - kavg) / kavg > 0.05:
        mid = (low + high) / 2
        G_mid,mid_mean = naive_config_gamma(mid,epsilon,n)
        if mid_mean < kavg:
            low = mid
        else:
            high = mid
    return G_mid,(low + high) / 2


def find_multi_k_binary_search(kavg,epsilon,n,net_type):
    high, low = n, kavg
    if net_type=='bet':
        high = n / (1+epsilon**2)
    mid = (low + high) / 2
    G_mean,mid_mean = naive_config_multi(mid,epsilon,n,net_type)
    while np.abs(mid_mean - kavg) / kavg > 0.05:
        mid = (low + high) / 2
        G_mid,mid_mean = naive_config_multi(mid,epsilon,n,net_type)
        if mid_mean < kavg:
            low = mid
        else:
            high = mid
    return G_mid,(low + high) / 2


def configuration_model_undirected_graph_mulit_type(kavg,epsilon,N,net_type,skewness):
    k_avg_graph,kavg_graph_max,kavg_graph_min = 0,2*kavg,kavg
    kavg_mid = (kavg_graph_max+kavg_graph_min)/2
    if N>100:
        while np.abs(kavg-k_avg_graph)/kavg>0.05:
            if net_type=='ig':
                wald_mu, wald_lambda = kavg, kavg / epsilon ** 2
                d = numpy.random.default_rng().wald(wald_mu,wald_lambda,N).astype(int)
            elif net_type=='bet':
                alpha_beta_dist,beta_beta_dist = (N - kavg * (1 + epsilon ** 2)) / (N * epsilon ** 2),((kavg - N) * (kavg - N + kavg * epsilon ** 2)) / (kavg * N * epsilon ** 2)
                d=(numpy.random.default_rng().beta(alpha_beta_dist,beta_beta_dist,N)*N).astype(int)
            elif net_type=='ln':
                mu_log_norm,sigma_log_norm = -(1 / 2) * np.log((1 + epsilon ** 2) / kavg ** 2),np.sqrt(2 * np.log(kavg) + np.log((1 + epsilon ** 2) / kavg ** 2))
                d = numpy.random.lognormal(mu_log_norm,sigma_log_norm,N).astype(int)
            elif net_type=='gam':
                theta, shape, k_avg_graph = epsilon ** 2 * kavg, 1 / epsilon ** 2, 0.0
                d = numpy.random.default_rng().gamma(shape, theta, N).astype(int)
            elif net_type=='bd':
                d1 = int((1 / 2)*(2 + skewness*epsilon - np.sqrt(4 + skewness ** 2)*epsilon)*kavg)
                d2 = int(kavg +  ((skewness + np.sqrt(4 + skewness**2))*epsilon*kavg)/2)
                fraction = (1 + skewness/np.sqrt(4 + skewness**2))/2
                G = random_bimodal_graph_skewed(d1,d2,N,fraction)
                return G,np.array([G.degree(n) for n in G.nodes()])
            elif net_type =='pl':
                a=epsilon
                if a < 5.0:
                    a_graph, b_graph = find_b_binary_search(kavg_mid, N, a)
                    d = configuration_model_powerlaw(a_graph, b_graph, N)
                else:
                    a_graph, b_graph = find_a_binary_search(k_avg_graph, N, a)
                    d = configuration_model_powerlaw(a_graph, b_graph, N)

            elif net_type == 'dg':  # New double Gaussian case
                # Generate degrees from double Gaussian distribution
                k1 = int((1 / 2)*(2 + skewness*epsilon - np.sqrt(4 + skewness ** 2)*epsilon)*kavg)
                k2 = int(kavg +  ((skewness + np.sqrt(4 + skewness**2))*epsilon*kavg)/2)
                fraction = (1 + skewness/np.sqrt(4 + skewness**2))/2
                sigma = 20
                d1 = np.random.normal(k1, sigma, int(N * fraction)).astype(int)
                d2 = np.random.normal(k2, sigma, N - len(d1)).astype(int)
                d = np.concatenate([d1, d2])


            # # Remove zeros from d
            # d = d[d != 0]
            # Replace zeros with ones
            d[d == 0] = 1
            if np.sum(d)%2!=0:
                d[int(len(d)*np.random.random())]+=1
            G = nx.configuration_model(d)
            G = nx.Graph(G)
            G.remove_edges_from(nx.selfloop_edges(G))
            graph_degrees= np.array([G.degree(n) for n in G.nodes()])
            k_avg_graph = np.mean(graph_degrees)
            if k_avg_graph<kavg:
                kavg_graph_min = k_avg_graph
            else:
                kavg_graph_max = k_avg_graph
            kavg_mid = (kavg_graph_min+kavg_graph_max)/2
            # Remove nodes with zero degree
            zero_degree_nodes = [n for n, d in G.degree() if d == 0]
            # G.remove_nodes_from(zero_degree_nodes)
        return G,np.array([G.degree(n) for n in G.nodes()])
    G,kavg_graph = find_multi_k_binary_search(kavg,epsilon,N,net_type)
    graph_degrees = np.array([G.degree(n) for n in G.nodes()])
    # # Remove nodes with zero degree
    # zero_degree_nodes = [n for n, d in G.degree() if d == 0]
    # G.remove_nodes_from(zero_degree_nodes)

    return G,graph_degrees


# def configuration_model_undirected_graph_mulit_type(kavg,epsilon,N,net_type,skewness):
#     eps_avg_graph,eps_graph_max,eps_graph_min = 0,2*epsilon,epsilon
#     eps_mid = (eps_graph_max+eps_graph_min)/2
#     k_avg_graph,kavg_graph_max,kavg_graph_min = 0,2*kavg,0.5*kavg
#     kavg_mid = (kavg_graph_max+kavg_graph_min)/2
#     sigma_max,sigma_min = kavg, 2
#     sigma_mid = (sigma_max+sigma_min)/2
#     median_graph,median_final = 0, kavg
#     if N>100:
#         while np.abs(median_graph-median_final)/median_final>0.05:
#
#             while np.abs(epsilon-eps_avg_graph)/epsilon>0.05:
#                 if net_type=='ig':
#                     wald_mu, wald_lambda = kavg, kavg / epsilon ** 2
#                     d = numpy.random.default_rng().wald(wald_mu,wald_lambda,N).astype(int)
#                 elif net_type=='bet':
#                     alpha_beta_dist,beta_beta_dist = (N - kavg * (1 + epsilon ** 2)) / (N * epsilon ** 2),((kavg - N) * (kavg - N + kavg * epsilon ** 2)) / (kavg * N * epsilon ** 2)
#                     d=(numpy.random.default_rng().beta(alpha_beta_dist,beta_beta_dist,N)*N).astype(int)
#                 elif net_type=='ln':
#                     mu_log_norm,sigma_log_norm = -(1 / 2) * np.log((1 + epsilon ** 2) / kavg ** 2),np.sqrt(2 * np.log(kavg) + np.log((1 + epsilon ** 2) / kavg ** 2))
#                     d = numpy.random.lognormal(mu_log_norm,sigma_log_norm,N).astype(int)
#                 elif net_type=='gam':
#                     theta, shape, k_avg_graph = epsilon ** 2 * kavg, 1 / epsilon ** 2, 0.0
#                     d = numpy.random.default_rng().gamma(shape, theta, N).astype(int)
#                 elif net_type=='bd':
#                     d1 = int((1 / 2)*(2 + skewness*epsilon - np.sqrt(4 + skewness ** 2)*epsilon)*kavg)
#                     d2 = int(kavg +  ((skewness + np.sqrt(4 + skewness**2))*epsilon*kavg)/2)
#                     fraction = (1 + skewness/np.sqrt(4 + skewness**2))/2
#                     G = random_bimodal_graph_skewed(d1,d2,N,fraction)
#                     return G,np.array([G.degree(n) for n in G.nodes()])
#                 elif net_type =='pl':
#                     a=epsilon
#                     if a < 5.0:
#                         a_graph, b_graph = find_b_binary_search(kavg_mid, N, a)
#                         d = configuration_model_powerlaw(a_graph, b_graph, N)
#                     else:
#                         a_graph, b_graph = find_a_binary_search(k_avg_graph, N, a)
#                         d = configuration_model_powerlaw(a_graph, b_graph, N)
#
#                 # elif net_type == 'dg':  # New double Gaussian case
#                 #     # Generate degrees from double Gaussian distribution
#                 #     k1 = int((1 / 2)*(2 + skewness*epsilon - np.sqrt(4 + skewness ** 2)*epsilon)*kavg)
#                 #     k2 = int(kavg +  ((skewness + np.sqrt(4 + skewness**2))*epsilon*kavg)/2)
#                 #     fraction = (1 + skewness/np.sqrt(4 + skewness**2))/2
#                 #     d1 = np.random.normal(k1, sigma, int(N * fraction)).astype(int)
#                 #     d2 = np.random.normal(k2, sigma, N - len(d1)).astype(int)
#                 #     d = np.concatenate([d1, d2])
#                 elif net_type == 'dg':  # New double Gaussian case
#                     # Set means and fraction for double Gaussian distribution
#                     k1 = int((1 / 2) * (2 + skewness * epsilon - np.sqrt(4 + skewness ** 2) * epsilon) * kavg_mid)
#                     k2 = int(kavg + ((skewness + np.sqrt(4 + skewness ** 2)) * epsilon * kavg_mid) / 2)
#                     fraction = (1 + skewness / np.sqrt(4 + skewness ** 2)) / 2
#
#                     # Generate degrees from first Gaussian
#                     d1 = np.random.normal(k1, sigma_mid, int(N * fraction)).astype(int)
#
#                     # Identify cutoff based on d1 to replicate in d2
#                     non_zero_d1 = d1[d1 > 0]
#                     if len(non_zero_d1) > 0:
#                         cutoff_d1 = min(non_zero_d1)
#                     else:
#                         cutoff_d1 = 1  # Fallback if all d1 values are zero or less, unlikely with large N
#
#                     # Generate degrees from second Gaussian with cutoff equivalent to d1's malformation
#                     d2 = np.random.normal(k2, sigma_mid, N - len(d1)).astype(int)
#                     d2 = d2[d2 >= k2 - (k1 - cutoff_d1)]  # Apply cutoff equivalent to the offset in d1
#
#                     # Ensure enough samples are present in d2, regenerate if needed
#                     while len(d2) < N - len(d1):
#                         additional_d2 = np.random.normal(k2, sigma_mid, (N - len(d1)) - len(d2)).astype(int)
#                         additional_d2 = additional_d2[additional_d2 >= k2 - (k1 - cutoff_d1)]
#                         d2 = np.concatenate([d2, additional_d2])
#
#                     # Combine the two Gaussian distributions
#                     d = np.concatenate([d1, d2])
#
#                 # # Remove zeros from d
#                 # d = d[d != 0]
#                 # Replace zeros with ones
#                 d[d == 0] = 1
#                 if np.sum(d)%2!=0:
#                     d[int(len(d)*np.random.random())]+=1
#                 G = nx.configuration_model(d)
#                 G = nx.Graph(G)
#                 G.remove_edges_from(nx.selfloop_edges(G))
#                 graph_degrees= np.array([G.degree(n) for n in G.nodes()])
#                 eps_avg_graph = np.std(graph_degrees)/np.mean(graph_degrees)
#                 if eps_avg_graph>epsilon:
#                     kavg_graph_min = np.mean(graph_degrees)
#                 else:
#                     kavg_graph_max = np.mean(graph_degrees)
#                 kavg_mid = (kavg_graph_max + kavg_graph_min) / 2
#                 # Remove nodes with zero degreel
#                 # zero_degree_nodes = [n for n, d in G.degree() if d == 0]
#                 # G.remove_nodes_from(zero_degree_nodes)
#             graph_degrees = np.array([G.degree(n) for n in G.nodes()])
#             median_graph = np.median(graph_degrees)
#             if median_graph < median_final:
#                 sigma_min = sigma_mid
#             else:
#                 sigma_max = sigma_mid
#             sigma_mid = (sigma_min + sigma_max) / 2
#             eps_avg_graph, eps_graph_max, eps_graph_min = 0, 2 * epsilon, epsilon
#             eps_mid = (eps_graph_max + eps_graph_min) / 2
#             k_avg_graph, kavg_graph_max, kavg_graph_min = 0, 2 * kavg, 0.5 * kavg
#             kavg_mid = (kavg_graph_max + kavg_graph_min) / 2
#         return G,np.array([G.degree(n) for n in G.nodes()])
#     G,kavg_graph = find_multi_k_binary_search(kavg,epsilon,N,net_type)
#     graph_degrees = np.array([G.degree(n) for n in G.nodes()])
#     # # Remove nodes with zero degree
#     # zero_degree_nodes = [n for n, d in G.degree() if d == 0]
#     # G.remove_nodes_from(zero_degree_nodes)
#
#     return G,graph_degrees


def configuration_model_undirected_graph_gamma(kavg,epsilon,N):
    if N>1000:
        theta,shape,k_avg_graph =epsilon**2*kavg,1/epsilon**2,0.0
        while np.abs(kavg-k_avg_graph)/kavg>0.05:
            d = numpy.random.default_rng().gamma(shape,theta,N).astype(int)
            if np.sum(d)%2!=0:
                d[int(len(d)*np.random.random())]+=1
            G = nx.configuration_model(d)
            G = nx.Graph(G)
            G.remove_edges_from(nx.selfloop_edges(G))
            k_avg_graph = np.mean([G.degree(n) for n in G.nodes()])
        return G
    G,kavg_graph = find_k_binary_search(kavg,epsilon,N)
    return G

# def configuration_model_dist_gen(avg,sigma,n):
#     class CustomDistribution(rv_discrete):
#         def _pmf(self, k, a,sigma):
#             return (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-(1/2)*((k-a)/sigma)**2)
#
#     custom_dist =CustomDistribution()
#     sequence = custom_dist.rvs(a=avg,sigma=sigma,size=n)
#     if (np.sum(sequence)%2!=0):
#         rand_bin = random.randint(1, n)
#     G = nx.configuration_model(sequence)
#     G = nx.Graph(G)
#     G.remove_edges_from(nx.selfloop_edges(G))
#     return G


def import_jason_network(file):
    a=[]
    with open(file,"r") as f:
        reader = csv.reader(f,delimiter=",")
        for i, line in enumerate(reader):
            a.append(line)
    return np.array(a)

def create_jason_edges(jason_edge):
    edges=[]
    for i in range(np.size(jason_edge,0)):
        x=list(filter(('').__ne__, jason_edge[i,:]))
        for j in x:
            if j is '':
                pass
            edges.append((i, int(j)))
    return edges

def jason_graph(file_name):
    csv_edges_jason=import_jason_network(file_name)
    edges=create_jason_edges(csv_edges_jason)
    G = nx.Graph()
    G.name = "myFile.csv"
    G.add_edges_from(edges)
    return G


if __name__ == '__main__':
    # N,k=400,100
    # interaction_strength = 100.5
    # d1_in, d1_out, d2_in, d2_out = 45,45,55,55
    # beta_inf,beta_sus = 0.1,0.1
    # susceptibility_avg = 1.0
    # infectability_avg = 1.0
    #
    # # G = nx.complete_graph(N)
    # G = nx.random_regular_graph(k, N)
    # susceptibility = 'bimodal'
    # infectability = 'bimodal'
    # choose_beta = lambda net_dist, avg, epsilon: np.random.normal(avg, epsilon * avg, N) \
    #     if net_dist == 'gauss' else np.random.gamma((avg / epsilon) ** 2, epsilon ** 2 / avg, N) \
    #     if net_dist == 'gamma' else np.zeros(N) if net_dist == 'z' else np.ones(
    #     N) if net_dist == 'ones' else netinithomo.bi_beta(N, epsilon, avg)
    #
    # epsilon_sus,epsilon_inf=0.1,0.1
    # beta_sus = choose_beta(susceptibility, susceptibility_avg, epsilon_sus)
    # beta_inf = choose_beta(infectability, infectability_avg, epsilon_inf)
    #
    # # random_bimodal_directed_graph(d1_in, d1_out, d2_in, d2_out, N)
    # G = netinithomo.intalize_lam_graph(G, N, beta_sus, beta_inf)
    # temp=nx.adjacency_data(G)
    # with open('json_data.json', 'w') as outfile:
    #     json.dump(temp, outfile)
    # export_network(G)
    # export_adj_network_cpp(G,'Gnull.csv')
    # print('This no love song')

    # fig=plt.figure()
    # G = rand_weighted_networks.bulid_weighted_network(N, k, 0.0, 0.0, interaction_strength)
    # G = netinithomo.set_graph_attriubute_DiGraph(G)
    # degree_in = [G.in_degree(n,'weight') for n in G.nodes()]
    # degree_out = [G.out_degree(n,'weight') for n in G.nodes()]
    # # plt.hist((degree_in, degree_out), bins=100)
    # plt.xlabel('Degree')
    # plt.ylabel('count')
    # plt.title('N='+str(N)+' k='+str(k)+' eta='+str(interaction_strength))
    # plt.hist(degree_out, bins=100)
    # fig.savefig('hist.png', dpi=400)
    # plt.show()

    # eta=np.linspace(0.1,2.5,5)
    # std_in,std_out=[],[]
    # for i in eta:
    #     G = rand_weighted_networks.bulid_weighted_network(N, k, 0.0, 0.0, i)
    #     G = netinithomo.set_graph_attriubute_DiGraph(G)
    #     degree_in = [G.in_degree(n, 'weight') for n in G.nodes()]
    #     degree_out = [G.out_degree(n, 'weight') for n in G.nodes()]
    #     std_in.append(np.std(degree_in))
    # plt.plot(eta,std_in)
    # plt.show()
    # G=random_bimodal_directed_graph(20, 50, 80, 50, 100, seed=None)
    # draw_basic_nx_g(G)
    # np.histogram(create_random_degree_sequence('gauss', 0.1, 108.5, 600),100)
    # plt.hist(create_random_degree_sequence('gauss', 0.1, 108.5, 600), bins='auto')
    # plt.show()
    # temp=configuration_model_directed_graph('gauss',0.1,0.1,108.5,600)
    # din,dout=create_random_degree_sequence('gauss', 0.1, 108.5, 600),create_random_degree_sequence('gauss', 0.1, 108.5, 600)
    # plt.hist((din,dout),bins='auto')
    # plt.show()
    # temp=equate_degree_sequence(din,dout,600)
    # plt.hist(temp, bins='auto')
    # plt.show()
    # degree=create_random_degree_sequence('gauss', 0.1, 108.5, 600)
    # temp=nx.utils.discrete_sequence(1000,degree)
    # v_in,v_out=[],[]
    # eps_in,eps_out=[0.1,0.1,0.1,0.1,0.1,0.1],[0.05,0.1,0.15,0.2,0.25,0.3]

    # class CustomDistribution(rv_discrete):
    #     def _pmf(self, k, a, b):
    #         return b * a / (1 + b * k) ** (a + 1)
    k,epsilon,N,net_type= 50,0.6,10000,'dg'
    skewness = 0.99
    a = epsilon
    # G = configuration_model_undirected_graph_gamma(k,epsilon,N)
    G,degree_sequence = configuration_model_undirected_graph_mulit_type(k,epsilon,N,net_type,skewness)
    plot_gamma_distribution(G,k,epsilon,N,net_type)
    print(np.median(degree_sequence))
    # custom_dist = CustomDistribution()
    # a,k,n=10.0,20,10000
    # b=1/(k*(a-1))
    # a_for_graph = 3.37
    # a_for_graph = 10.0
    # a_graph, b_graph = find_a_binary_search(k, n, a)
    # k,n=10,100000
    # a_array,b_array=[0.5,1.1],[]
    # G_array,k_array=[],[]
    # for a in a_array:
    #     a_graph, b_graph = find_b_binary_search(k, n, a)
    #     G= configuration_model_powerlaw(a_graph, b_graph, n)
    #     k_avg_graph = np.mean([G.degree(n) for n in G.nodes()])
    #     while (np.abs(k_avg_graph-k)/k>0.05):
    #         if a<5.0:
    #             a_graph, b_graph = find_b_binary_search(k, n, a)
    #         else:
    #             a_graph, b_graph = find_a_binary_search(k, n, a)
    #         G= configuration_model_powerlaw(a_graph, b_graph, n)
    #         k_avg_graph = np.mean([G.degree(n) for n in G.nodes()])
    #     b_array.append(b_graph)
    #     k_avg_graph = np.mean([G.degree(n) for n in G.nodes()])
    #     k_array.append(k_avg_graph)
    #     G_array.append(G)
    # plot_configuration_model_powerlaw_multi(G_array,a_array,b_array,n,custom_dist,k_array)
    # plot_configuration_model_powerlaw(G,a,b,n,custom_dist,k_avg_graph)
    # degree_sequence = sorted((d for n, d in G.degree()), reverse=True)
    # plt.hist(degree_sequence, bins=np.arange(1, 50)-0.5, density=True, alpha=0.6, color='g', label='Histogram of Generated degree')
    # plt.title('Power-law Degree Distribution')
    # plt.xlabel('Degree (k)')
    # plt.ylabel('Probability')
    # x = np.arange(1, 50)
    # plt.plot(x, custom_dist.pmf(x,a=a, b=b), 'ro-', label='PMF', markersize=8)
    # plt.show()
    # eps_in,eps_out=[1.0],[1.0]
    # for ein,eout in zip(eps_in,eps_out):
    #     G=configuration_model_directed_graph('gauss_c', ein, eout, 50, 400)
    #     degree_in=[G.in_degree(n) for n in G.nodes()]
    #     degree_out=[G.out_degree(n) for n in G.nodes()]
    #     # v_in.append(np.std(degree_in)/np.mean(degree_in))
    #     # v_out.append(np.std(degree_out)/np.mean(degree_out))
    #     # plt.hist((degree_in,degree_out),bins=500)
    #     plt.hist(degree_in, bins=100)
    #     # plt.hist((degree_out),bins=400)
    #     plt.show()
    # avg_degree,epsilon,N=[20,20,20],[0.9786772558354417,0.5,0.2],[10000,10000,10000]
    # for eps,avg,n in zip(epsilon,avg_degree,N):
    #     # G=configuration_model_undirected_graph(eps,avg,n)
    #     G = configuration_model_undirected_graph_exp(avg, n)
    #     degree=np.array([G.degree(n) for n in G.nodes()])
    #     # v_in.append(np.std(degree_in)/np.mean(degree_in))
    #     # v_out.append(np.std(degree_out)/np.mean(degree_out))
    #     # plt.hist((degree_in,degree_out),bins=500)
    #     plt.xlabel('k')
    #     plt.ylabel('p(k)')
    #     plt.title(r'Histogram Gauss for N={} $\sigma={}$, <k>={}'.format(n,round(np.std(degree),2),round(np.mean(degree),2)))
    #     plt.hist(degree, bins=200,density=True)
    #     plt.savefig('hist_degree_N{}_eps_{}_kavg_{}.png'.format(n,round(eps,2),avg),dpi=500)
    #     plt.show()
    # print(*v_in, sep = ", ")
    # print(*v_out, sep = ", ")
    # degrees_in = [G.in_degree(n) for n in G.nodes()]
    # degrees_out = [G.out_degree(n) for n in G.nodes()]
    # plt.hist((degrees_in, degrees_out), bins='auto')
    # plt.show()
    # temp=import_jason_network('jason_trans_file_no_degree.csv')
    # temp_edge=create_jason_edges(temp)
    # G=jason_graph('jason_trans_file_no_degree.csv')
    # print('This no love song')
# k=4
# d1=78
# d2=66
# N=796

#This set of commands is used in order to create a random homogenous graph and draw it

# G=nx.random_regular_graph(k,N)
# draw_basic_nx_g(G)
# export_network(G)


#This set of commands is used in order to create a bimodal graph and draw it

# G=random_bimodal_graph(d1 , d2, N, seed=None)
# export_network(G)
# draw_basic_nx_g(G)
