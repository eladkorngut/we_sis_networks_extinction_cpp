#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <random>
#include <sstream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <tuple>
#include <numeric>
#include <experimental/filesystem>
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "queue"
#include "sstream"
#include <unistd.h>



// Declaration of functions
void remove_infected_node(int node,std::vector<int>& infected_node,std::vector<int> &positions);
void add_susceptible_node(int node,int m,std::vector<std::vector<int>>& susceptible_node,std::vector<int>& positions);
void decrement_susc_neighbs(int k_in,int k_out,std::vector<int>& neighbs_in,std::vector<int>& neighbs_out,
                            std::vector<std::vector<int>>& susceptible_nodes,std::vector<int>& infected_neighbors_in,
                            std::vector<int>& infected_neighbors_out,std::vector<int>& s_m,std::vector<int>& positions,std::vector<int>& sigma);
//void increment_susc_neighbs(int k_in,int k_out,std::vector<int>& neighbs_in,std::vector<int>& neighbs_out,std::vector<std::vector<int>>& susceptible_nodes,
//                            std::vector<int>& infected_neighbs_int,std::vector<int>& infected_neighbs_out,std::vector<int>& s_m,std::vector<int>& positions,std::vector<int>& sigma);
void remove_susceptible_node(int node,int k,std::vector<std::vector<int>>& susceptible_nodes,std::vector<int>& positions);
void add_infected_node(int node,std::vector<int>& infected_node,std::vector<int>& positions);
void read_in_neighborslist(int N,std::string& filename,std::vector<std::vector<int>>& Adjlist,std::vector<int>& degree,std::fstream& Adjfile);
void read_parameters(std::string& filename,std::fstream& parameters_file,int N,int sims,int it,int k,double x,double lam,double Alpha,double Beta,int network_number,
                     double mf_solution,int new_trajectory_bin,int k_max);


// Start of auxiliary functions

int read_single_int_parameter(std::stringstream& ss,std::string& word){
    getline(ss, word, ',');
    ss.clear();
    return std::stoi(word);

}

double read_single_double_parameter(std::stringstream& ss,std::string& word) {
    getline(ss, word, ',');
    ss.clear();
    return std::stof(word);
}

std::string read_single_string_parameter(std::stringstream& ss,std::string& word) {
    getline(ss, word, ',');
    ss.clear();
    return word;
}

std::tuple<int,int,int,int,double,double,double,double,double,int,double,double,int,std::string & ,std::string & >
        read_parameters(std::string& filename,std::fstream& parameters_file) {
    // Read network parameters from python files
    std::string line;
    std::stringstream ss;
    std::string word;
    int N,sims,it,k,jump,network_number,new_trajectory_bin;
    double x,lam,Alpha,Beta,tau,mf_solution;
    std::string prog,dir_path;
    parameters_file.open("/home/elad/we_sis_networks_extinction_cpp/cparameters.txt");
    while (getline(parameters_file, line)) {
        ss.clear();
        ss.str(line);
        N = read_single_int_parameter(ss,word);
        sims = read_single_int_parameter(ss,word);
        it = read_single_int_parameter(ss,word);
        k = read_single_int_parameter(ss,word) ;
        x = read_single_double_parameter(ss,word);
        lam = read_single_double_parameter(ss,word);
        jump = read_single_int_parameter(ss,word);
        Alpha = read_single_double_parameter(ss,word);
        Beta = read_single_double_parameter(ss,word);
        network_number = read_single_int_parameter(ss,word);
        tau = read_single_double_parameter(ss,word);
        mf_solution = read_single_double_parameter(ss,word);
        new_trajectory_bin = read_single_int_parameter(ss,word);
        prog=read_single_string_parameter(ss,word);
        dir_path=read_single_string_parameter(ss,word);
    }
    return std::tuple<int,int,int,int,double,double,int,double,double,int,double,double,int,std::string & ,std::string&>
            (N,sims,it,k,x,lam,jump,Alpha,Beta,network_number,tau,mf_solution,new_trajectory_bin,prog,dir_path);
}

void remove_infected_node(int node,std::vector<int> &infected_node,std::vector<int> &positions){
    // Remove the node from the list by swapping the node with the last node in the list and removing the last node in the list. This will not affect the positioning of the others
    infected_node[positions[node]] = infected_node.back();
    positions[infected_node.back()] = positions[node];
    infected_node.pop_back();
}

void add_susceptible_node(int node,int m,std::vector<std::vector<int>>& susceptible_node,std::vector<int>& positions){
    positions[node] = susceptible_node[m].size();
    susceptible_node[m].push_back(node);
}


void decrement_susc_neighbs(int k_in,int k_out,std::vector<int> &neighbs_in,std::vector<int> &neighbs_out,
                            std::vector<std::vector<int>> &net_susceptible_nodes,std::vector<int> &infected_neighbors_in,
                            std::vector<int> &infected_neighbors_out,std::vector<int> &s_m,std::vector<int> &positions
                            ,std::vector<int> &sigma,int node,std::vector<int> &inf_neigh,int temp_step_count,int temp_t_count){
    int neighb,m2;
    for (int j=0;j<k_out;j++){
        // Change neighbours from s_m class to s_m{m-1} class. need to be for the ingoing neighbours
        neighb = neighbs_out[j];
        infected_neighbors_in[neighbs_out[j]]--;
        if (sigma[neighb]==0){
            m2 = infected_neighbors_in[neighb];
            remove_susceptible_node(neighb,m2+1,net_susceptible_nodes,positions);
            s_m[m2+1]--;
            add_susceptible_node(neighb,m2,net_susceptible_nodes,positions);
            s_m[m2]++;
        }

    }
    for (int j=0;j<k_in;j++){
        infected_neighbors_out[neighbs_in[j]]--; // Each neighbor of recovered node has one less infected neighbor
    }
}

//void decrement_susc_neighbs(int k_in,int k_out,std::vector<int>& neighbs_in,std::vector<int>& neighbs_out,
//                            std::vector<std::vector<int>>& net_susceptible_nodes,std::vector<int>& infected_neighbors_in,
//                            std::vector<int>& infected_neighbors_out,std::vector<int>& s_m,std::vector<int>& positions
//        ,std::vector<int>& sigma,int node,std::vector<int> &inf_neigh){
//    int neighb,m2;
//    for (int j=0;j<k_out;j++){
//        // Each neighbor of recovered node has one less infected neighbor
//        infected_neighbors_in[neighbs_out[j]]--;
//    }
//    for (int j=0;j<k_in;j++){
//        // Change neighbours from s_m class to s_m{m-1} class. need to be for the ingoing neighbours
//        neighb = neighbs_in[j];
//        infected_neighbors_out[neighb]--;
//        if (sigma[neighb]==0){
//            m2 = infected_neighbors_in[neighb];
//            remove_susceptible_node(neighb,m2+1,net_susceptible_nodes,positions);
//            s_m[m2+1]--;
//            add_susceptible_node(neighb,m2,net_susceptible_nodes,positions);
//            s_m[m2]++;
//        }
//    }
//}

void remove_susceptible_node(int node,int m, std::vector<std::vector<int>> &net_susceptible_nodes,std::vector<int> &positions){

    // Remove the node from the list by swapping the node with the last node in the list. This will not affect the
    // positioning of the other nodes once the position of the swapped last node is accounted for
    int temp,temp2;
    temp = net_susceptible_nodes[m].back();
    temp2 = net_susceptible_nodes[m][positions[node]];
//    net_susceptible_nodes[m][positions[node]] = temp;
    net_susceptible_nodes[m][positions[node]] = net_susceptible_nodes[m].back();  // Move last node in list to nodes position
    positions[net_susceptible_nodes[m].back()] = positions[node]; // account for this change in the position vector
    net_susceptible_nodes[m].pop_back(); // remove the last node from the list
}


//void increment_susc_neighbs(int k_in,int k_out,std::vector<int>& neighbs_in,std::vector<int>& neighbs_out,
//                            std::vector<std::vector<int>> &net_susceptible_nodes,std::vector<int>& infected_neighbs_in,
//                            std::vector<int>& infected_neighbs_out,std::vector<int>& s_m,std::vector<int>& positions,
//                            std::vector<int>& sigma){
//    int neighb,m;
//    for (int i = 0; i < k_out; i++) {
//        // Each neighbour of infected node has one more infected neighbour
//        infected_neighbs_out[neighbs_out[i]]++;
//    }
//    for (int i=0; i<k_in;i++){
//        neighb = neighbs_in[i];
//        infected_neighbs_out[neighb]++;
//        if (sigma[neighb]==0){
//            m = infected_neighbs_in[neighb];
//            remove_susceptible_node(neighb,m-1,net_susceptible_nodes,positions);
//            s_m[m-1]--;
//            add_susceptible_node(neighb,m,net_susceptible_nodes,positions);
//            s_m[m]++;
//        }
//    }
//}

void increment_susc_neighbs(int k_in,int k_out,std::vector<int> &neighbs_in,std::vector<int> &neighbs_out,
                            std::vector<std::vector<int>> &net_suceptible_nodes,std::vector<int> &infected_neighbors_in,
                            std::vector<int> &infected_neighbors_out,std::vector<int> &s_m,std::vector<int> &positions,
                            std::vector<int> &sigma){
    int neighb,m2;
    for (int j=0;j<k_in;j++){
        infected_neighbors_out[neighbs_in[j]]++;
    }
    for (int j=0;j<k_out;j++){
        neighb = neighbs_out[j];
        infected_neighbors_in[neighb]++;
        if (sigma[neighb]==0){
            m2 = infected_neighbors_in[neighb];
            remove_susceptible_node(neighb,m2-1,net_suceptible_nodes,positions);
            s_m[m2-1]--;
            add_susceptible_node(neighb,m2,net_suceptible_nodes,positions);
            s_m[m2]++;
        }
    }
}


void add_infected_node(int node,std::vector<int> &infected_nodes,std::vector<int> &positions){

    positions[node]= infected_nodes.size();
    infected_nodes.push_back(node);
//    positions[node]= infected_nodes.size()-1;
}


void read_in_neighborslist(int N,std::string& filename,std::vector<std::vector<int>> &Adjlist,std::vector<int> &degrees,std::fstream &Adjfile,int &k_max) {
    // Read network data from python files

    Adjfile.open("/home/elad/we_sis_networks_extinction_cpp/"+filename);
    std::string line;
    std::stringstream ss;
    int count_degree = 0;
    std::string word;
    if (Adjfile.is_open()) {
        while (std::getline(Adjfile, line)) {
            ss.clear();
            ss.str(line);
            std::getline(ss, word, ',');
            degrees[count_degree] = std::stoi(word);
            if (degrees[count_degree]>k_max)
                k_max = degrees[count_degree];
            while (getline(ss, word, ',')) {
                Adjlist[count_degree].push_back(std::stoi(word));

            }
            count_degree++;
        }
    }
    Adjfile.close();
}


//void increment_neighbours(int node,int k,std::vector<int>& neighbs,std::vector<int>& infected_neighbs_in){
//    for (int i=0;i<k;i++)
//        infected_neighbs_in[neighbs[i]]++;
//}

void increment_neighbours(int k,std::vector<int> &neighbs,std::vector<int> &infected_neighbs_in){
    for (int i=0;i<k;i++)
        infected_neighbs_in[neighbs[i]]++;
}

void intalized_network_random_infection(int N,int inital_infecteons,std::list< std::vector<int>>& sigma,std::list<std::vector<int>>& infected_node
                                        ,std::list<std::vector<int>>&  positions){
    int count(0),node;
    for (std::list< std::vector<int>>::iterator it_sigma=sigma.begin(),it_infected_node=infected_node.begin(),
                 it_positions=positions.begin();it_sigma!=sigma.end();++it_sigma,++it_infected_node,++it_positions){
        count=0;
        while (inital_infecteons>count){
            node = rand()%N;
            if ((*it_sigma)[node]==0)
            {
                count++;
                (*it_sigma)[node]=1;
                (*it_positions)[node] = it_infected_node->size();
                it_infected_node->push_back(node);
            }
        }
    }
}

//void intalize_infected_neighbours(int N,std::list< std::vector<int>>& sigma,std::list< std::vector<int>>& infected_neighbors_in,
//                                  std::list< std::vector<int>>& infected_neighbors_out,std::vector<int>& degrees_in
//                                  ,std::vector<int>& degrees_out,std::vector<std::vector<int>>& Adjlist_in,std::vector<std::vector<int>>& Adjlist_out){
//    for (std::list< std::vector<int>>::iterator it_infected_neighbors_in=infected_neighbors_in.begin(),
//                 it_infected_neighbors_out=infected_neighbors_out.begin(),it_sigma=sigma.begin();it_infected_neighbors_in!=infected_neighbors_in.end()
//            ;++it_infected_neighbors_in,++it_infected_neighbors_out,++it_sigma){
//        for (int i=0;i<N;i++){
//            // if node is in state 1, increment neighbours
//            if ((*it_sigma)[i]==1){
//                increment_neighbours(i,degrees_out[i],Adjlist_out[i],*it_infected_neighbors_in);
//                increment_neighbours(i,degrees_in[i],Adjlist_in[i],*it_infected_neighbors_out);
//            }
//        }
//    }
//}


void intalize_infected_neighbours(int N,std::list< std::vector<int>> &infected,std::list< std::vector<int>> &infected_neighbors_in,
                                  std::list< std::vector<int>> &infected_neighbors_out,std::vector<int> &degrees_in
        ,std::vector<int>& degrees_out,std::vector<std::vector<int>>& Adjlist_in,std::vector<std::vector<int>>& Adjlist_out){
    for (std::list< std::vector<int>>::iterator it_infected_neighbors_in=infected_neighbors_in.begin(),
                 it_infected_neighbors_out=infected_neighbors_out.begin(),it_infected=infected.begin();it_infected_neighbors_in!=infected_neighbors_in.end()
            ;++it_infected_neighbors_in,++it_infected_neighbors_out,++it_infected){
        for (int i=0;i<it_infected->size();i++){
                increment_neighbours(degrees_out[(*it_infected)[i]],Adjlist_out[(*it_infected)[i]],*it_infected_neighbors_in);
                increment_neighbours(degrees_in[(*it_infected)[i]],Adjlist_in[(*it_infected)[i]],*it_infected_neighbors_out);
        }
    }
}


void intalized_sus(int N,std::list< std::vector<std::vector<int>>> &susceptible_nodes,std::list<std::vector<int>> &positions,std::list< std::vector<int>> &sigma,
                   std::list<std::vector<int>> &s_m,std::list< std::vector<int>> &infected_neighbors_in){
    std::list< std::vector<std::vector<int>>>::iterator it_susceptible_nodes=susceptible_nodes.begin();
    std::list<std::vector<int>>::iterator it_positions=positions.begin();
    std::list< std::vector<int>>::iterator it_sigma=sigma.begin();
    std::list<std::vector<int>>::iterator it_s_m=s_m.begin();
    std::list< std::vector<int>>::iterator it_infected_neighbors_in=infected_neighbors_in.begin();
    while (it_susceptible_nodes!=susceptible_nodes.end()){
        for (int i=0;i<N;i++){
            if ((*it_sigma)[i]==0){
                (*it_positions)[i] = (*it_susceptible_nodes)[(*it_infected_neighbors_in)[i]].size();
                (*it_susceptible_nodes)[(*it_infected_neighbors_in)[i]].push_back(i);
                (*it_s_m)[(*it_infected_neighbors_in)[i]]++;
            }
        }
        ++it_susceptible_nodes;
        ++it_positions;
        ++it_sigma;
        ++it_s_m;
        ++it_infected_neighbors_in;
    }

}

//void intalize_SI_connections(int N,std::list<int>& SI,std::list< std::vector<int>>& sigma,std::list< std::vector<int>>& infected_neighbors_in,
//                             std::list< std::vector<int>>& infected_neighbors_out,std::vector<int>& degrees_out){
//    std::list<int>::iterator it_SI=SI.begin();
//    auto it_sigma=sigma.begin();
//    auto it_infected_neighbors_in=infected_neighbors_in.begin();
//    auto it_infected_neighbors_out=infected_neighbors_out.begin();
//    while (it_sigma!=sigma.end()){
//        for (int i=0;i<N;i++){
//            if ((*it_sigma)[i]==0){
//                (*it_SI)+=(*it_infected_neighbors_in)[i];
//            }
//            else{
//                (*it_SI)+=degrees_out[i]-(*it_infected_neighbors_out)[i];
//            }
//        }
//        ++it_SI;
//        ++it_sigma;
//        ++it_infected_neighbors_in;
//    }
//}
void  intalize_SI_connections(std::list<int>& SI,std::list<std::vector<int>>& s_m,int k_max){
    std::list<std::vector<int>>::iterator it_s_m=s_m.begin();
    std::list<int>::iterator it_SI=SI.begin();
    while(it_SI!=SI.end()){
        for (int i=0;i<=k_max;i++) {
            (*it_SI) += (*it_s_m)[i] * i;
        }
        ++it_SI;
        ++it_s_m;
    }
}

void inital_networks_stat(int N,double x,int sims,std::list<int>& num_inf,std::list<std::vector<int>>& infected_node,
                          std::list< std::vector<int>>& infected_neighbors_in,std::list< std::vector<int>>& infected_neighbors_out,
                          std::list< std::vector<int>>& sigma,std::vector<int>& degrees_in,std::vector<int>& degrees_out,
                          std::vector<std::vector<int>>& Adjlist_in,std::vector<std::vector<int>>& Adjlist_out,std::list<std::vector<int>>&  positions
                          ,std::list< std::vector<std::vector<int>>>& susceptible_nodes,std::list<int>& SI,std::list<std::vector<int>>& s_m,int k_max){
    int inital_infecteons(int(x*N));
    num_inf = std::list<int>(sims,inital_infecteons);
    susceptible_nodes = std::list< std::vector<std::vector<int>>>(sims,std::vector<std::vector<int>>(k_max+1,std::vector<int>(0)));
    intalized_network_random_infection(N,inital_infecteons,sigma,infected_node,positions);
//    std::mt19937 generator(std::random_device{}()); //random number generator
//    intalize_infected_neighbours(N,sigma,infected_neighbors_in,infected_neighbors_out,degrees_in,degrees_out,Adjlist_in,Adjlist_out);
    intalize_infected_neighbours(N,infected_node,infected_neighbors_in,infected_neighbors_out,degrees_in,degrees_out,Adjlist_in,Adjlist_out);
    intalized_sus(N,susceptible_nodes,positions,sigma,s_m,infected_neighbors_in);
//    intalize_SI_connections(N,SI,sigma,infected_neighbors_in,infected_neighbors_out,degrees_out);
    intalize_SI_connections(SI,s_m,k_max);
}


void GillespieMC(int steps_c,std::list<int> &num_inf,std::list<int>::iterator &it_num_inf,
                 std::list<double>::iterator &it_avec_sum,std::list<double> &avec_sum,
                 std::list<double>::iterator &it_t,std::list<double> &t,
                 std::list<std::vector<int>> &infected_node,std::list<std::vector<int>>::iterator &it_infected_node,
                 std::list< std::vector<int>> &infected_neighbors_in,std::list< std::vector<int>>::iterator &it_infected_neighbors_in,
                 std::list< std::vector<int>> &infected_neighbors_out,std::list< std::vector<int>>::iterator &it_infected_neighbors_out,
                 std::list<std::vector<int>> &sigma,std::list<std::vector<int>>::iterator &it_sigma,
                 std::list<std::vector<int>> &s_m,std::list<std::vector<int>>::iterator &it_s_m,
                 std::list<std::vector<int>> &positions,std::list<std::vector<int>>::iterator &it_positions,
                 std::list< std::vector<std::vector<int>>> &susceptible_nodes,std::list< std::vector<std::vector<int>>>::iterator& it_susceptible_nodes,
                 std::list<int> &SI,std::list<int>::iterator &it_SI,std::list<bool> &skip,std::list<bool>::iterator &it_skip,
                 double Alpha,double tau,int k_max,double death,double Beta,std::vector<int> &degrees_in,std::vector<int> &degrees_out,
                 std::vector<std::vector<int>> &Adjlist_in,std::vector<std::vector<int>> &Adjlist_out,std::list<double>::iterator &it_weights){
    int temp_count_steps(0),temp_count_num_inf(0);
    int m_in,m_out,node;
    double r1,r2,s_m_sum;

    while (0<steps_c){
        //The initialized the iterators of the different list
        temp_count_steps++;
        temp_count_num_inf=0;

        it_num_inf=num_inf.begin();
        it_avec_sum=avec_sum.begin();
        it_t=t.begin();
        it_infected_node=infected_node.begin();
        it_infected_neighbors_in=infected_neighbors_in.begin();
        it_infected_neighbors_out=infected_neighbors_out.begin();
        it_sigma=sigma.begin();
        it_s_m=s_m.begin();
        it_positions=positions.begin();
        it_susceptible_nodes =susceptible_nodes.begin();
        it_SI = SI.begin();
        it_skip = skip.begin();


        // For all networks draw the exponentially distributed time step
        while (it_num_inf!=num_inf.end()){
            if (*it_skip){
                ++it_num_inf;
                ++it_avec_sum;
                ++it_t;
                ++it_infected_node;
                ++it_infected_neighbors_in;
                ++it_infected_neighbors_out;
                ++it_sigma;
                ++it_s_m;
                ++it_positions;
                ++it_susceptible_nodes;
                ++it_SI;
                ++it_skip;
                continue;
            }
            temp_count_num_inf++;
            r1 = rand()/(double(RAND_MAX));
            r2 = std::log(1.0/(double(rand())/double(RAND_MAX)));
            *it_avec_sum = (*it_num_inf)*Alpha+(*it_SI)*Beta;
            *it_t = *it_t + r2/(*it_avec_sum);
            if (*it_t>tau){
                (*it_skip)=true;
                steps_c =steps_c-1;
            }
//            *it_t =*(std::prev(it_t)) + r2/(*it_avec_sum);

            //Pick a node, change its state and update the transition rates
            if (r1<(*it_num_inf)*Alpha/(*it_avec_sum)) {
                // Recover an infected node
                node = (*it_infected_node)[int(rand())%(*it_num_inf)]; //randomly choose infected node
//                m = (*it_infected_neighbors)[node]; //number of infected neighbors
                m_in = (*it_infected_neighbors_in)[node]; //number of infected neighbors that infecting the node
                m_out = (*it_infected_neighbors_out)[node]; //number of infected neighbors that the node infect
                (*it_sigma)[node]=0;
                (*it_num_inf)--;
                //Remove node from infected list and add to susceptible list
                remove_infected_node(node,*it_infected_node,*it_positions);
                if ((*it_num_inf)==0){
                    //Advance to the next simulation in the list
                    ++it_num_inf;
                    ++it_avec_sum;
                    ++it_t;
                    ++it_infected_node;
                    ++it_infected_neighbors_in;
                    ++it_infected_neighbors_out;
                    ++it_sigma;
                    ++it_s_m;
                    ++it_positions;
                    ++it_susceptible_nodes;
                    ++it_SI;
                    continue;
                }
//                (*it_SI) = 2*m-degrees[node];
                (*it_SI) = (*it_SI) + m_in + m_out - degrees_out[node];

//                remove_infected_node(node,*it_infected_node,*it_positions);
//                add_susceptible_node(node,m,*it_susceptible_nodes,*it_positions);
//                if ((*it_num_inf)<=0){
//                    continue;
//                }
                add_susceptible_node(node,m_in,*it_susceptible_nodes,*it_positions);


                // Increment number of susceptible nodes with m infected nodes
                (*it_s_m)[m_in]++;


                decrement_susc_neighbs(degrees_in[node],degrees_out[node],Adjlist_in[node],Adjlist_out[node],*it_susceptible_nodes,
                                       *it_infected_neighbors_in,*it_infected_neighbors_out,*it_s_m,*it_positions,*it_sigma,node,*it_infected_node,temp_count_steps,temp_count_num_inf);


            }
            else{
                //Infected a susceptible node
                // Pick an m class with probability proportional to the total number of SI links in that m class i.e. s_m[m]*m
                r1= rand()/(double(RAND_MAX));
                r1=r1*(*it_SI);
                s_m_sum = 0;
//                for (int m1=0;m1<=k_max;m1++){
                for (int m1=1;m1<=k_max;m1++){
                    s_m_sum = s_m_sum + (*it_s_m)[m1]*m1;
                    if (r1<=s_m_sum){
                        // choose a node with m infected neighbours at random
                        node = (*it_susceptible_nodes)[m1][rand()%(*it_s_m)[m1]];
                        break;
                    }
                }
//                m = (*it_infected_neighbors)[node];
                m_in = (*it_infected_neighbors_in)[node];
                m_out = (*it_infected_neighbors_out)[node];
                (*it_sigma)[node] = 1;
                *it_num_inf = (*it_num_inf)+1;
//                *it_SI = *it_SI - 2*m - degrees[node];
//                (*it_SI) = 2*m_in+2*m_out-degrees_in[node]-degrees_out[node];


                // Remove the node from susceptible list and add to infected list
                remove_susceptible_node(node,m_in,*it_susceptible_nodes,*it_positions);
                add_infected_node(node,*it_infected_node,*it_positions);

                //Decrement number of susceptible nodes
                (*it_s_m)[m_in]--;

                // change the neighbours to other lists and adjust s_m

//                increment_susc_neighbs(degrees_in[node],degrees_out[node],Adjlist_in[node],Adjlist_out[node],
//                                       *it_susceptible_nodes,*it_infected_neighbors_in,*it_infected_neighbors_out,
//                                       *it_s_m,*it_positions,*it_sigma);
//                if (susceptible_nodes.size()<198){
//                    bool stop= true;
//                }
                increment_susc_neighbs(degrees_in[node],degrees_out[node],Adjlist_in[node],Adjlist_out[node],
                                       *it_susceptible_nodes,*it_infected_neighbors_in,*it_infected_neighbors_out,
                                       *it_s_m,*it_positions,*it_sigma);

                (*it_SI) = (*it_SI) + degrees_out[node] - m_in - m_out;
            }
            //Advance to the next simulation in the list
            ++it_num_inf;
            ++it_avec_sum;
            ++it_t;
            ++it_infected_node;
            ++it_infected_neighbors_in;
            ++it_infected_neighbors_out;
            ++it_sigma;
            ++it_s_m;
            ++it_positions;
            ++it_susceptible_nodes;
            ++it_SI;
            ++it_skip;
        } //end one time step

        temp_count_num_inf=0;

        it_num_inf=num_inf.begin();
        it_avec_sum=avec_sum.begin();
        it_t=t.begin();
        it_infected_node=infected_node.begin();
        it_infected_neighbors_in=infected_neighbors_in.begin();
        it_infected_neighbors_out=infected_neighbors_out.begin();
        it_sigma=sigma.begin();
        it_s_m=s_m.begin();
        it_positions=positions.begin();
        it_susceptible_nodes =susceptible_nodes.begin();
        it_SI = SI.begin();

        while (it_num_inf!=num_inf.end()){
            temp_count_num_inf++;
            if (*it_num_inf<=0){
                death = death + (*it_weights);
                it_num_inf=num_inf.erase(it_num_inf);
                it_avec_sum=avec_sum.erase(it_avec_sum);
                it_t=t.erase(it_t);
                it_infected_node=infected_node.erase(it_infected_node);
                it_infected_neighbors_in=infected_neighbors_in.erase(it_infected_neighbors_in);
                it_infected_neighbors_out=infected_neighbors_out.erase(it_infected_neighbors_out);
                it_sigma=sigma.erase(it_sigma);
                it_s_m=s_m.erase(it_s_m);
                it_positions=positions.erase(it_positions);
                it_susceptible_nodes=susceptible_nodes.erase(it_susceptible_nodes);
                it_SI=SI.erase(it_SI);
                steps_c =steps_c-1;
                if (it_num_inf==num_inf.end()){
                    break;
                }
                continue;
            }
            ++it_num_inf;
            ++it_avec_sum;
            ++it_t;
            ++it_infected_node;
            ++it_infected_neighbors_in;
            ++it_infected_neighbors_out;
            ++it_sigma;
            ++it_s_m;
            ++it_positions;
            ++it_susceptible_nodes;
            ++it_SI;
            ++it_skip;
        }
        std::cout <<"Time "<<  temp_count_steps <<" "<< "Simulation "<<temp_count_num_inf<<std::endl;
        temp_count_num_inf=0;
    }
}


int main() {
    int m_in,m_out,node;  //k_max maximal degree, N number of nodes, sims number of simulations
    double death(0),r1,r2,s_m_sum;
    std::fstream Adjfile_in,Adjfile_out,parametersfile;
    std::string filename_in("Adjin_0.txt"), filename_out("Adjout_0.txt"), parametersname("parameters.csv");
    auto parameter_list=read_parameters(parametersname,parametersfile);
    int N=std::get<0>(parameter_list),sims=std::get<1>(parameter_list),it=std::get<2>(parameter_list),
                k=std::get<3>(parameter_list),jump=std::get<6>(parameter_list),network_number=std::get<9>(parameter_list),new_trajectory_bin=std::get<12>(parameter_list);
    double inital_inf_percent=std::get<4>(parameter_list),lam=std::get<5>(parameter_list),Alpha=std::get<7>(parameter_list)
                ,Beta=std::get<8>(parameter_list),tau=std::get<10>(parameter_list),mf_solution=std::get<11>(parameter_list);
                std::string dir_path(std::get<13>(parameter_list)),prog=std::get<14>(parameter_list);
    int steps_c(sims),k_max(0),k_max_out(0); // Need to find ouf which dimension i np.size(n,1) refers to
    std::list<double> weights(sims,1.0/double(sims)),wg(sims,1.0/double(sims)),avec_sum(sims,0),t(sims,0);
    std::list<double>::iterator it_weights(weights.begin()),it_wg(wg.begin()),it_avec_sum(avec_sum.begin());
    std::vector<int> degrees_in(N,0),degrees_out(N,0);
    std::vector<std::vector<int>> Adjlist_in(N,std::vector<int>()),Adjlist_out(N,std::vector<int>());
    std::list<std::vector<int>> sigma(sims,std::vector<int>(N,0)),infected_node(sims,std::vector<int>()),
    ng(steps_c,std::vector<int>(N,0)),positions(sims,std::vector<int>(N,0));
    std::list<std::vector<int>>::iterator it_sigma(sigma.begin()),it_infected_node(infected_node.begin()),it_ng(ng.begin()),it_positions(positions.begin());
    std::list< std::vector<std::vector<int>>> susceptible_nodes(steps_c,std::vector<std::vector<int>>());
    std::list< std::vector<std::vector<int>>>::iterator it_susceptible_nodes(susceptible_nodes.begin());
    std::list<int> num_inf(sims,0),SI(sims,0);
    std::list<int>::iterator it_num_inf(num_inf.begin()),it_SI(SI.begin());
    std::list<double>::iterator it_t(t.begin());
    std::list< std::vector<int>> infected_neighbors_in(steps_c,std::vector<int>(N,0)),infected_neighbors_out(steps_c,std::vector<int>(N,0));
    std::list< std::vector<int>>::iterator it_infected_neighbors_in(infected_neighbors_in.begin()),it_infected_neighbors_out(infected_neighbors_out.begin());
    std::list<bool> skip(sims, false);
    std::list<bool>::iterator it_skip=skip.begin();


    //Read the python's network structure, what are the different connections

    read_in_neighborslist(N,filename_in,Adjlist_in,degrees_in, Adjfile_in,k_max);
    read_in_neighborslist(N,filename_out,Adjlist_out,degrees_out,Adjfile_out,k_max_out);

    // Intailze the different networks list such as the number of infected nodes and so on
    std::list<std::vector<int>> s_m(steps_c,std::vector<int>(k_max+1,0));
    std::list<std::vector<int>>::iterator it_s_m(s_m.begin());
    inital_networks_stat(N,inital_inf_percent,sims,num_inf,infected_node,infected_neighbors_in,infected_neighbors_out,sigma,
                         degrees_in,degrees_out,Adjlist_in,Adjlist_out,positions,susceptible_nodes,SI,s_m,k_max);
    int temp_count_steps(0),temp_count_num_inf(0);
    
    GillespieMC(steps_c,num_inf,it_num_inf,it_avec_sum,avec_sum,it_t,t,infected_node,it_infected_node,infected_neighbors_in,it_infected_neighbors_in,
            infected_neighbors_out,it_infected_neighbors_out,sigma,it_sigma,s_m,it_s_m,positions,it_positions,susceptible_nodes,it_susceptible_nodes,
            SI,it_SI,skip,it_skip,Alpha,tau,k_max,death,Beta,degrees_in,degrees_out,Adjlist_in,Adjlist_out,it_weights);

//    while (0<steps_c){
//        //The initialized the iterators of the different list
//        temp_count_steps++;
//        temp_count_num_inf=0;
//
//        it_num_inf=num_inf.begin();
//        it_avec_sum=avec_sum.begin();
//        it_t=t.begin();
//        it_infected_node=infected_node.begin();
//        it_infected_neighbors_in=infected_neighbors_in.begin();
//        it_infected_neighbors_out=infected_neighbors_out.begin();
//        it_sigma=sigma.begin();
//        it_s_m=s_m.begin();
//        it_positions=positions.begin();
//        it_susceptible_nodes =susceptible_nodes.begin();
//        it_SI = SI.begin();
//
//
//        // For all networks draw the exponentially distributed time step
//        while (it_num_inf!=num_inf.end()){
//            temp_count_num_inf++;
//            r1 = rand()/(double(RAND_MAX));
//            r2 = std::log(1.0/(double(rand())/double(RAND_MAX)));
//            *it_avec_sum = (*it_num_inf)*Alpha+(*it_SI)*Beta;
//            *it_t = *it_t + r2/(*it_avec_sum);
////            *it_t =*(std::prev(it_t)) + r2/(*it_avec_sum);
//
//            //Pick a node, change its state and update the transition rates
//            if (r1<(*it_num_inf)*Alpha/(*it_avec_sum)) {
//                // Recover an infected node
//                node = (*it_infected_node)[int(rand())%(*it_num_inf)]; //randomly choose infected node
////                m = (*it_infected_neighbors)[node]; //number of infected neighbors
//                m_in = (*it_infected_neighbors_in)[node]; //number of infected neighbors that infecting the node
//                m_out = (*it_infected_neighbors_out)[node]; //number of infected neighbors that the node infect
//                (*it_sigma)[node]=0;
//                (*it_num_inf)--;
//                //Remove node from infected list and add to susceptible list
//                remove_infected_node(node,*it_infected_node,*it_positions);
//                if ((*it_num_inf)==0){
//                    //Advance to the next simulation in the list
//                    ++it_num_inf;
//                    ++it_avec_sum;
//                    ++it_t;
//                    ++it_infected_node;
//                    ++it_infected_neighbors_in;
//                    ++it_infected_neighbors_out;
//                    ++it_sigma;
//                    ++it_s_m;
//                    ++it_positions;
//                    ++it_susceptible_nodes;
//                    ++it_SI;
//                    continue;
//                }
////                (*it_SI) = 2*m-degrees[node];
//                (*it_SI) = (*it_SI) + m_in + m_out - degrees_out[node];
//
////                remove_infected_node(node,*it_infected_node,*it_positions);
////                add_susceptible_node(node,m,*it_susceptible_nodes,*it_positions);
////                if ((*it_num_inf)<=0){
////                    continue;
////                }
//                add_susceptible_node(node,m_in,*it_susceptible_nodes,*it_positions);
//
//
//                // Increment number of susceptible nodes with m infected nodes
//                (*it_s_m)[m_in]++;
//
//
//                decrement_susc_neighbs(degrees_in[node],degrees_out[node],Adjlist_in[node],Adjlist_out[node],*it_susceptible_nodes,
//                                       *it_infected_neighbors_in,*it_infected_neighbors_out,*it_s_m,*it_positions,*it_sigma,node,*it_infected_node,temp_count_steps,temp_count_num_inf);
//
//
//            }
//            else{
//                //Infected a susceptible node
//                // Pick an m class with probability proportional to the total number of SI links in that m class i.e. s_m[m]*m
//                r1= rand()/(double(RAND_MAX));
//                r1=r1*(*it_SI);
//                s_m_sum = 0;
////                for (int m1=0;m1<=k_max;m1++){
//                for (int m1=1;m1<=k_max;m1++){
//                    s_m_sum = s_m_sum + (*it_s_m)[m1]*m1;
//                    if (r1<=s_m_sum){
//                        // choose a node with m infected neighbours at random
//                        node = (*it_susceptible_nodes)[m1][rand()%(*it_s_m)[m1]];
//                        break;
//                    }
//                }
////                m = (*it_infected_neighbors)[node];
//                m_in = (*it_infected_neighbors_in)[node];
//                m_out = (*it_infected_neighbors_out)[node];
//                (*it_sigma)[node] = 1;
//                *it_num_inf = (*it_num_inf)+1;
////                *it_SI = *it_SI - 2*m - degrees[node];
////                (*it_SI) = 2*m_in+2*m_out-degrees_in[node]-degrees_out[node];
//
//
//                // Remove the node from susceptible list and add to infected list
//                remove_susceptible_node(node,m_in,*it_susceptible_nodes,*it_positions);
//                add_infected_node(node,*it_infected_node,*it_positions);
//
//                //Decrement number of susceptible nodes
//                (*it_s_m)[m_in]--;
//
//                // change the neighbours to other lists and adjust s_m
//
////                increment_susc_neighbs(degrees_in[node],degrees_out[node],Adjlist_in[node],Adjlist_out[node],
////                                       *it_susceptible_nodes,*it_infected_neighbors_in,*it_infected_neighbors_out,
////                                       *it_s_m,*it_positions,*it_sigma);
////                if (susceptible_nodes.size()<198){
////                    bool stop= true;
////                }
//                increment_susc_neighbs(degrees_in[node],degrees_out[node],Adjlist_in[node],Adjlist_out[node],
//                                       *it_susceptible_nodes,*it_infected_neighbors_in,*it_infected_neighbors_out,
//                                       *it_s_m,*it_positions,*it_sigma);
//
//                (*it_SI) = (*it_SI) + degrees_out[node] - m_in - m_out;
//            }
//            //Advance to the next simulation in the list
//            ++it_num_inf;
//            ++it_avec_sum;
//            ++it_t;
//            ++it_infected_node;
//            ++it_infected_neighbors_in;
//            ++it_infected_neighbors_out;
//            ++it_sigma;
//            ++it_s_m;
//            ++it_positions;
//            ++it_susceptible_nodes;
//            ++it_SI;
//        } //end one time step
//
//        temp_count_num_inf=0;
//
//        it_num_inf=num_inf.begin();
//        it_avec_sum=avec_sum.begin();
//        it_t=t.begin();
//        it_infected_node=infected_node.begin();
//        it_infected_neighbors_in=infected_neighbors_in.begin();
//        it_infected_neighbors_out=infected_neighbors_out.begin();
//        it_sigma=sigma.begin();
//        it_s_m=s_m.begin();
//        it_positions=positions.begin();
//        it_susceptible_nodes =susceptible_nodes.begin();
//        it_SI = SI.begin();
//
//        while (it_num_inf!=num_inf.end()){
//            temp_count_num_inf++;
//            if (*it_num_inf<=0){
//                 death = death + (*it_weights);
//                it_num_inf=num_inf.erase(it_num_inf);
//                it_avec_sum=avec_sum.erase(it_avec_sum);
//                it_t=t.erase(it_t);
//                it_infected_node=infected_node.erase(it_infected_node);
//                it_infected_neighbors_in=infected_neighbors_in.erase(it_infected_neighbors_in);
//                it_infected_neighbors_out=infected_neighbors_out.erase(it_infected_neighbors_out);
//                it_sigma=sigma.erase(it_sigma);
//                it_s_m=s_m.erase(it_s_m);
//                it_positions=positions.erase(it_positions);
//                it_susceptible_nodes=susceptible_nodes.erase(it_susceptible_nodes);
//                it_SI=SI.erase(it_SI);
//                steps_c =steps_c-1;
//                if (it_num_inf==num_inf.end()){
//                    break;
//                }
//                continue;
//            }
//            ++it_num_inf;
//            ++it_avec_sum;
//            ++it_t;
//            ++it_infected_node;
//            ++it_infected_neighbors_in;
//            ++it_infected_neighbors_out;
//            ++it_sigma;
//            ++it_s_m;
//            ++it_positions;
//            ++it_susceptible_nodes;
//            ++it_SI;
//        }
//        std::cout <<"Time "<<  temp_count_steps <<" "<< "Simulation "<<temp_count_num_inf<<std::endl;
//        temp_count_num_inf=0;
//    }
    std::cout << "Hello, World!" << std::endl;
    return 0;
}