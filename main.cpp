#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <random>

void remove_infected_node(int node,std::vector<int>& infected_node,std::vector<int> positions){}
void add_susceptible_node(int node,int m,std::vector<std::vector<int>>& susceptible_node,std::vector<int>& positions){}
void decrement_susc_neighbs(int node,int k,std::vector<int>& neighbors,std::vector<std::vector<int>>& susceptible_nodes,std::vector<int>& infected_neighbors,
                            std::vector<int>& s_m,std::vector<int>& positions,std::vector<int>& sigma){}
void remove_susceptible_node(int node,int k,std::vector<std::vector<int>>& susceptible_nodes,std::vector<int>& positions){}
void add_infected_node(int node,std::vector<int>& infected_node,std::vector<int>& positions){}


int main() {
    int k_max;  //probably the maximal degree
    int sims =5;
    int N;
    std:std::vector< std::vector<int> > n(sims,std::vector<int>(0,N)); // Initial infection state row simulations and column are nodes
    std::list<double> weights;
    std::list<double>::iterator it_weights(weights.begin());
    std::vector<int> degrees;
    float tau(1.0);
    std::vector<std::vector<int>> Adjlist(N,std::vector<int>(N,1));
    std::list< std::vector<int>> sigma;
    std::list<std::vector<int>>::iterator it_sigma(sigma.begin());
    std::list<std::vector<int>> infected_node;
    std::list<std::vector<int>>::iterator it_infected_node(infected_node.begin());
    std::list< std::vector<std::vector<int>>> susceptible_nodes;
    std::list< std::vector<std::vector<int>>>::iterator it_susceptible_nodes(susceptible_nodes.begin());
    int steps_c(n.size()); // Need to find ouf which dimension i np.size(n,1) refers to
//    int sims((*n.begin()).size());
    // Because I later on remove some elements it might be best to use list instead of vector to define the variables below
    std::list<std::vector<int>> ng(steps_c,std::vector<int>(N,0));
    std::list<std::vector<int>>::iterator it_ng(ng.begin());
    std::list<std::vector<int>> s_m(steps_c,std::vector<int>(k_max+1,0)); //Number of infected neighbours of each node
    std::list<std::vector<int>>::iterator it_s_m(s_m.begin());
    std::list<float> wg(steps_c,0);
    std::list<float>::iterator it_wg(wg.begin());
    std::list<float> avec_sum(sims,0);
    std::list<float>::iterator it_avec_sum(avec_sum.begin());
    double death(0);
    float r1(0),r2(0);
    float Alpha(1.0);
    std::list<int> num_inf(sims,0);
    std::list<int>::iterator it_num_inf(num_inf.begin());
    float beta(0);
    std::list<int> SI_links(sims,0);
    std::list<int>::iterator it_SI_links(SI_links.begin());
    std::list<float> t(sims,0);
    std::list<float>::iterator it_t(t.begin());
    std::list< std::vector<int>> infected_neighbors(steps_c,std::vector<int>(N,0));
    std::list< std::vector<int>>::iterator it_infected_neighbors(infected_neighbors.begin());
    std::list<std::vector<int>>  positions(sims,std::vector<int>(N,0));
    std::list<std::vector<int>>::iterator it_positions(positions.begin());
    int node(0);
    std::list<int> SI(sims,0);
    std::list<int>::iterator it_SI(SI.begin());
    int m;
    float s_m_sum;
    while (1){
        //The initialized the iterators of the different list



        // For all networks draw the exponentially distributed time step
        for (int i=0; i<sims;i++){
            r1 = float(rand())/float(RAND_MAX);
            r2 = std::log(float(rand())/float(RAND_MAX));
            avec_sum[i] = num_inf[i]*Alpha+SI_links[i]*beta;
            t[i] =t[i]+ r2/avec_sum[i];

            //Pick a node, change its state and update the transition rates
            if (r1<num_inf[i]*Alpha/avec_sum[i]) {
                // Recover an infected node
                node = infected_node[i][rand()%num_inf[i]]; //randomly choose infected node
                m = infected_neighbors[i][node]; //number of infected neighbors
                sigma[i][node]=0;
                num_inf[i]--;
                SI_links[i] = 2*m-degrees[node];

                //Remove node from infected list and add to susceptible list
                remove_infected_node(node,infected_node[i],positions[i]);
                add_susceptible_node(node,m,susceptible_nodes[i],positions[i]);

                // Increment number of susceptible nodes with m infected nodes
                s_m[i][m]++;

                decrement_susc_neighbs(node,degrees[node],Adjlist[node],susceptible_nodes[i],infected_neighbors[i],s_m[i],positions[i],sigma[i]);


            }
            else{
                //Infected a susceptible node

                // Pick an m class with probability proprtional to the total number of SI links in that m class i.e. s_m[m]*m

                r1= float(rand())/float(RAND_MAX);
                r1=r1*SI[i];
                s_m_sum = 0;
                for (int m1=0;m1<=k_max;m1++){
                    s_m_sum = s_m_sum + s_m[i][m1]*m1;
                    if (r1<=s_m_sum){
                        // choose a node with m infected neighbours at random
                        node = susceptible_nodes[i][m1][rand()%s_m[i][m1]];
                        break;
                    }
                }
                m = infected_neighbors[i][node];
                sigma[i][node] = 1;
                num_inf[i] = num_inf[i]+1;
                SI[i] = SI[i] - 2*m - degrees[node];

                // Remove the node from susceptible list and add to infected list
                remove_susceptible_node(node,m,susceptible_nodes[i],positions[i]);
                add_infected_node(node,infected_node[i],positions[i]);

            }
        } //end one time step
        for (int i=0;i<num_inf.size();i++){
            if (num_inf[i]<=0){
                death = death + weights[i];

            }
        }


    }





    std::cout << "Hello, World!" << std::endl;
    return 0;
}
