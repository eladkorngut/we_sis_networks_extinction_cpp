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


// Declaration of functions
void remove_infected_node(int node,std::vector<int>& infected_node,std::vector<int> positions);
void add_susceptible_node(int node,int m,std::vector<std::vector<int>>& susceptible_node,std::vector<int>& positions);
void decrement_susc_neighbs(int node,int k,std::vector<int>& neighbors,std::vector<std::vector<int>>& susceptible_nodes,std::vector<int>& infected_neighbors,
                            std::vector<int>& s_m,std::vector<int>& positions,std::vector<int>& sigma);
void increment_susc_neighbs(int node,int k,std::vector<int>& neighbs,std::vector<std::vector<int>>& susceptible_nodes,
                            std::vector<int>& infected_neighbs,std::vector<int>& s_m,std::vector<int>& positions,std::vector<int>& sigma);
void remove_susceptible_node(int node,int k,std::vector<std::vector<int>>& susceptible_nodes,std::vector<int>& positions);
void add_infected_node(int node,std::vector<int>& infected_node,std::vector<int>& positions);
void read_in_neighborslist(int N,int k_max,std::string filename,std::vector<std::vector<int>>& Adjlist,std::fstream& Adjfile);
void read_parameters(int N,int sims,int it,int k,float x,float lam,float Alpha,float Beta,int network_number,std::string& infile,
                     float mf_solution,int new_trajectory_bin);

void remove_infected_node(int node,std::vector<int>& infected_node,std::vector<int> positions){
    infected_node[positions[node]] = infected_node.back();
    positions[infected_node.back()] = positions[node];
    infected_node.pop_back();
}

void add_susceptible_node(int node,int m,std::vector<std::vector<int>>& susceptible_node,std::vector<int>& positions){
    positions[node] = susceptible_node[m].size();
    susceptible_node[m].push_back(node);
}
void decrement_susc_neighbs(int node,int k,std::vector<int>& neighbs,std::vector<std::vector<int>>& susceptible_nodes,std::vector<int>& infected_neighbors,
                            std::vector<int>& s_m,std::vector<int>& positions,std::vector<int>& sigma){
    int neighb,m2;
    for (int j=0;j<k;j++){

        neighb = neighbs[j];

        // Each neighbor of recovered node has one less infected neighbor
        infected_neighbors[neighb]--;

        // Change neighbourz from s_m class to s_m{m-1} class
         if (sigma[neighb]==0){
             m2 = infected_neighbors[neighb];
             remove_susceptible_node(neighb,m2,susceptible_nodes,positions);
             s_m[m2+1]--;
             add_susceptible_node(neighb,m2,susceptible_nodes,positions);
             s_m[m2]++;
         }
    }
}


void increment_susc_neighbs(int node,int k,std::vector<int>& neighbs,std::vector<std::vector<int>>& susceptible_nodes,
                            std::vector<int>& infected_neighbs,std::vector<int>& s_m,std::vector<int>& positions,std::vector<int>& sigma){
    int neighb,m;
    for (int (i) = 0; (i) < k; ++(i)) {

        neighb = neighbs[i];

        // Each neighbour of infected node has one more infected neighbour
        infected_neighbs[neighb]++;

        if (sigma[neighb]==0){
            m = infected_neighbs[neighb];
            remove_susceptible_node(neighb,m-1,susceptible_nodes,positions);
            s_m[m-1]--;
            add_susceptible_node(neighb,m,susceptible_nodes,positions);
            s_m[m]++;
        }
    }
}


void remove_susceptible_node(int node,int m, std::vector<std::vector<int>>& susceptible_nodes,std::vector<int>& positions){

    // Remove the node from the list by swapping the node with the last node in the list. This will not affect the
    // positioning of the other nodes once the position of the swapped last node is accounted for

    susceptible_nodes[m][positions[node]] = susceptible_nodes[m].back();  // Move last node in list to nodes position
    positions[susceptible_nodes[m].back()] = positions[node]; // account for this change in the position vector
    susceptible_nodes[m].pop_back(); // remove the last node from the list
}


void add_infected_node(int node,std::vector<int>& infected_nodes,std::vector<int>& positions){

    positions[node]= infected_nodes.size();
    infected_nodes.push_back(node);
}


void read_in_neighborslist(int N,std::string& filename,std::vector<std::vector<int>>& Adjlist,std::vector<int>& degrees,std::fstream& Adjfile){
    // Read network data from python files
    Adjfile.open(filename);
    for (int i=0;i<N;i++){
        Adjfile >> degrees[i];
        for (int j=0;j<degrees[i];j++){
            Adjfile>> Adjlist[i][j];
        }
    }

}


void increment_neighbours(int node,int k,std::vector<int>& neighbs,std::vector<int>& infected_neighbs){
    for (int i=0;i<k;i++)
        infected_neighbs[neighbs[i]]++;
}


void inital_networks_stat(int N,float x,std::string& filename,std::fstream& Adjfile,int sims,std::list<int>& num_inf,std::list<std::vector<int>>& infected_node,
                               std::list< std::vector<int>>& infected_neighbors,std::list< std::vector<int>>& sigma,std::vector<int>& degrees
                               ,std::vector<std::vector<int>>& Adjlist,std::list<std::vector<int>>&  positions,std::list< std::vector<std::vector<int>>>& susceptible_nodes,
                                std::list<int>& SI,std::list<std::vector<int>>& s_m){
    int inital_infecteons(int(x*N));
    num_inf = std::list<int>(int(x*N),0);
    sigma = std::list<std::vector<int>>(sims,std::vector<int>(N,0));
    susceptible_nodes = std::list< std::vector<std::vector<int>>>(sims,std::vector<std::vector<int>>(N,std::vector<int>(0)));
//    std::mt19937 generator(std::random_device{}()); //random number generator

    int count(0),node;
    for (std::list< std::vector<int>>::iterator it_sigma=sigma.begin(),it_infected_node=infected_node.begin(),
            it_positions=positions.begin();it_sigma!=sigma.end();++it_sigma,++it_infected_node,++it_positions){
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

    for (std::list< std::vector<int>>::iterator it_infected_neighbors=infected_neighbors.begin(), it_sigma=sigma.begin();it_infected_neighbors!=infected_neighbors.end();++it_infected_neighbors,++it_sigma){
        for (int i=0;i<N;i++){
            // if node is in state 1, increment neighbours
            if ((*it_sigma)[i]==1){
                increment_neighbours(i,degrees[i],Adjlist[i],*it_infected_neighbors);
            }
        }
    }

    count = 0;
    std::list< std::vector<std::vector<int>>>::iterator it_susceptible_nodes=susceptible_nodes.begin();
    std::list<std::vector<int>>::iterator it_positions=positions.begin();
    std::list< std::vector<int>>::iterator it_sigma=sigma.begin();
    std::list<std::vector<int>>::iterator it_s_m=s_m.begin();
    std::list< std::vector<int>>::iterator it_infected_neighbors=infected_neighbors.begin();
    while (it_susceptible_nodes!=susceptible_nodes.end()){
        for (int i=0;i<N;i++){
            if ((*it_sigma)[i]==0){
                (*it_positions)[i] = (*it_susceptible_nodes)[(*it_infected_neighbors)[i]].size();
                (*it_susceptible_nodes)[(*it_infected_neighbors)[i]].push_back(i);
                (*it_s_m)[(*it_infected_neighbors)[i]]++;
            }
        }
        ++it_susceptible_nodes;
        ++it_positions;
        ++it_sigma;
        ++it_s_m;
        ++it_infected_neighbors;
    }


    SI = std::list<int>(sims,0);
    std::list<int>::iterator it_SI=SI.begin();
    it_sigma=sigma.begin();
    it_infected_neighbors=infected_neighbors.begin();
    while (it_sigma!=sigma.end()){
        for (int i=0;i<N;i++){
            if ((*it_sigma)[i]==0)
            {
                (*it_SI)+=(*it_infected_neighbors)[i];
            }
        }
        ++it_SI;
        ++it_sigma;
        ++it_infected_neighbors;
    }
}



int main() {
    int k_max;  //probably the maximal degree
    int sims =5;
    int N;
//    std:std::vector< std::vector<int> > n(sims,std::vector<int>(0,N)); // Initial infection state row simulations and column are nodes
    std::list<double> weights;
    std::list<double>::iterator it_weights(weights.begin());
    std::vector<int> degrees;
    double tau(1.0);
    std::vector<std::vector<int>> Adjlist(N,std::vector<int>(N,1));
    std::list< std::vector<int>> sigma;
    std::list<std::vector<int>>::iterator it_sigma(sigma.begin());
    std::list<std::vector<int>> infected_node;
    std::list<std::vector<int>>::iterator it_infected_node(infected_node.begin());
    std::list< std::vector<std::vector<int>>> susceptible_nodes;
    std::list< std::vector<std::vector<int>>>::iterator it_susceptible_nodes(susceptible_nodes.begin());
    int steps_c(sims); // Need to find ouf which dimension i np.size(n,1) refers to
//    int sims((*n.begin()).size());
    // Because I later on remove some elements it might be best to use list instead of vector to define the variables below
    std::list<std::vector<int>> ng(steps_c,std::vector<int>(N,0));
    std::list<std::vector<int>>::iterator it_ng(ng.begin());
    std::list<std::vector<int>> s_m(steps_c,std::vector<int>(k_max+1,0)); //Number of infected neighbours of each node
    std::list<std::vector<int>>::iterator it_s_m(s_m.begin());
    std::list<double> wg(steps_c,0);
    std::list<double>::iterator it_wg(wg.begin());
    std::list<double> avec_sum(sims,0);
    std::list<double>::iterator it_avec_sum(avec_sum.begin());
    double death(0);
    double r1(0),r2(0);
    double Alpha(1.0);
    std::list<int> num_inf(sims,0);
    std::list<int>::iterator it_num_inf(num_inf.begin());
    double beta(0);
    std::list<double> t(sims,0);
    std::list<double>::iterator it_t(t.begin());
    std::list< std::vector<int>> infected_neighbors(steps_c,std::vector<int>(N,0));
    std::list< std::vector<int>>::iterator it_infected_neighbors(infected_neighbors.begin());
    std::list<std::vector<int>>  positions(sims,std::vector<int>(N,0));
    std::list<std::vector<int>>::iterator it_positions(positions.begin());
    int node(0);
    std::list<int> SI(sims,0);
    std::list<int>::iterator it_SI(SI.begin());
    int m;
    double s_m_sum;
    while (steps_c<=0){
        //The initialized the iterators of the different list

        it_num_inf=num_inf.begin();
        it_avec_sum=avec_sum.begin();
        it_t=t.begin();
        it_infected_node=infected_node.begin();
        it_infected_neighbors=infected_neighbors.begin();
        it_sigma=sigma.begin();
        it_s_m=s_m.begin();
        it_positions=positions.begin();
        it_susceptible_nodes =susceptible_nodes.begin();
        it_SI = SI.begin();


        // For all networks draw the exponentially distributed time step
        while (it_num_inf!=num_inf.end()){
            r1 = double(rand())/double(RAND_MAX);
            r2 = std::log(double(rand())/double(RAND_MAX));
            *it_avec_sum = (*it_num_inf)*Alpha+(*it_SI)*beta;
            *it_t =*it_t+ r2/(*it_avec_sum);

            //Pick a node, change its state and update the transition rates
            if (r1<(*it_num_inf)*Alpha/(*it_avec_sum)) {
                // Recover an infected node
                node = (*it_infected_node)[rand()%(*it_num_inf)]; //randomly choose infected node
                m = (*it_infected_neighbors)[node]; //number of infected neighbors
                (*it_sigma)[node]=0;
                (it_num_inf)--;
                (*it_SI) = 2*m-degrees[node];

                //Remove node from infected list and add to susceptible list
                remove_infected_node(node,*it_infected_node,*it_positions);
                add_susceptible_node(node,m,*it_susceptible_nodes,*it_positions);

                // Increment number of susceptible nodes with m infected nodes
                (*it_s_m)[m]++;

                decrement_susc_neighbs(node,degrees[node],Adjlist[node],*it_susceptible_nodes,*it_infected_neighbors,*it_s_m,*it_positions,*it_sigma);


            }
            else{
                //Infected a susceptible node

                // Pick an m class with probability proprtional to the total number of SI links in that m class i.e. s_m[m]*m

                r1= double(rand())/double(RAND_MAX);
                r1=r1*(*it_SI);
                s_m_sum = 0;
                for (int m1=0;m1<=k_max;m1++){
                    s_m_sum = s_m_sum + (*it_s_m)[m1]*m1;
                    if (r1<=s_m_sum){
                        // choose a node with m infected neighbours at random
                        node = (*it_susceptible_nodes)[m1][rand()%(*it_s_m)[m1]];
                        break;
                    }
                }
                m = (*it_infected_neighbors)[node];
                (*it_sigma)[node] = 1;
                *it_num_inf = (*it_num_inf)+1;
                *it_SI = *it_SI - 2*m - degrees[node];

                // Remove the node from susceptible list and add to infected list
                remove_susceptible_node(node,m,*it_susceptible_nodes,*it_positions);
                add_infected_node(node,*it_infected_node,*it_positions);

                //Decrement number of susceptible nodes
                (*it_s_m)[m]--;

                // change the neighbours to other lists and adjust s_m

                increment_susc_neighbs(node,degrees[node],Adjlist[node],*it_susceptible_nodes,*it_infected_neighbors,*it_s_m,*it_positions,*it_sigma);

            }
            //Advance to the next simulation in the list
            ++it_num_inf;
            ++it_avec_sum;
            ++it_t;
            ++it_infected_node;
            ++it_infected_neighbors;
            ++it_sigma;
            ++it_s_m;
            ++it_positions;
            ++it_susceptible_nodes;
            ++it_SI;
        } //end one time step

        it_num_inf=num_inf.begin();
        it_avec_sum=avec_sum.begin();
        it_t=t.begin();
        it_infected_node=infected_node.begin();
        it_infected_neighbors=infected_neighbors.begin();
        it_sigma=sigma.begin();
        it_s_m=s_m.begin();
        it_positions=positions.begin();
        it_susceptible_nodes =susceptible_nodes.begin();
        it_SI = SI.begin();


        while (it_num_inf!=num_inf.end()){
            if (*it_num_inf<=0){
                death = death + (*it_weights);
                num_inf.erase(it_num_inf);
                avec_sum.erase(it_avec_sum);
                t.erase(it_t);
                infected_node.erase(it_infected_node);
                infected_neighbors.erase(it_infected_neighbors);
                sigma.erase(it_sigma);
                s_m.erase(it_s_m);
                positions.erase(it_positions);
                susceptible_nodes.erase(it_susceptible_nodes);
                SI.erase(it_SI);
                steps_c =steps_c-1;
            }
        }


//        for (int i=0;i<num_inf.size();i++){
//            if (*it_num_inf<=0){
//                death = death + *it_weights;
//
//
//            }
//        }


    }





    std::cout << "Hello, World!" << std::endl;
    return 0;
}
