#include <iostream>
#include <fstream>
#include <vector>
#include <lemon/smart_graph.h>

using namespace std;
using namespace lemon;

struct Server {
    vector<vector<unsigned int>> consumption_rate;
    vector<unsigned int> transition_costs;
    unsigned int amount_servers;

    Server(vector<vector<unsigned int>> cr, vector<unsigned int> tc, unsigned int amount) {
        consumption_rate = move(cr);
        transition_costs = move(tc);
        amount_servers = amount;
    }
};

void generate_graph(SmartDigraph &g, SmartDigraph::NodeMap<long> &imbalances1, SmartDigraph::NodeMap<long> &imbalances2,
                    SmartDigraph::ArcMap<unsigned int> &cost, SmartDigraph::ArcMap<unsigned int> &capacity,
                    SmartDigraph::ArcMap<unsigned int> &capacity1, SmartDigraph::ArcMap<unsigned int> &capacity2,
                    vector<Server> &servers, vector<unsigned int> &demands) {
    unsigned int amount_servers = 0;
    unsigned int d_k = 0;
    for(auto it_server = servers.begin(); it_server != servers.end(); ++it_server) {
        amount_servers += it_server->amount_servers;
        d_k += it_server->amount_servers * (it_server->transition_costs.size() - 1);
    }

    vector<SmartDigraph::Node> sources;
    vector<SmartDigraph::Node> sinks;

    SmartDigraph::Node a_0 = g.addNode();
    imbalances1.set(a_0, amount_servers);
    sources.push_back(a_0);

    SmartDigraph::Node b_0 = g.addNode();
    imbalances1.set(b_0, -1 * amount_servers);
    sinks.push_back(b_0);

    for(auto it_demands = demands.begin(); it_demands != demands.end(); ++it_demands) {
        SmartDigraph::Node a_k = g.addNode();
        imbalances2.set(a_k, d_k + *it_demands);
        sources.push_back(a_k);

        SmartDigraph::Node b_k = g.addNode();
        imbalances2.set(b_k, -1 * (d_k + *it_demands));
        sinks.push_back(b_k);
    }

    SmartDigraph::Node old_u;
    vector<SmartDigraph::Node> old_l;
    for(auto it_servers = servers.begin(); it_servers != servers.end(); ++it_servers) {
        for(int j = 0; j != demands.size(); ++j) {
            SmartDigraph::Node u_k = g.addNode();
            vector<SmartDigraph::Node> current_l;
            if(j != 0) {
                //Connect this timestep with prior one
                SmartDigraph::Arc edge = g.addArc(old_u, u_k);
                capacity.set(edge, it_servers->amount_servers);
                capacity1.set(edge, it_servers->amount_servers);
                capacity2.set(edge, 0);
                cost.set(edge, it_servers->consumption_rate[0][j-1]);
            }
            old_u = u_k;
            for(int k = 0; k != it_servers->transition_costs.size(); ++k) {
                SmartDigraph::Node l_km = g.addNode();
                SmartDigraph::Node l_kma = g.addNode();
                current_l.push_back(l_kma);
                if(j == 0 && k == it_servers->transition_costs.size() - 1) {
                    SmartDigraph::Arc edge = g.addArc(sources[0], l_km);
                    capacity.set(edge, it_servers->amount_servers);
                    capacity1.set(edge, it_servers->amount_servers);
                    capacity2.set(edge, 0);
                    cost.set(edge, 0);

                    edge = g.addArc(l_km, u_k);
                    capacity.set(edge, it_servers->amount_servers);
                    capacity1.set(edge, it_servers->amount_servers);
                    capacity2.set(edge, 0);
                    cost.set(edge, it_servers->transition_costs[k]);
                }
                if(j != 0) {
                    SmartDigraph::Arc edge = g.addArc(l_km, u_k);
                    capacity.set(edge, it_servers->amount_servers);
                    capacity1.set(edge, it_servers->amount_servers);
                    capacity2.set(edge, 0);
                    cost.set(edge, it_servers->transition_costs[k]);

                    edge = g.addArc(u_k, l_km);
                    capacity.set(edge, it_servers->amount_servers);
                    capacity1.set(edge, it_servers->amount_servers);
                    capacity2.set(edge, 0);
                    cost.set(edge, 0);

                    edge = g.addArc(l_km, sinks[j]);
                    capacity.set(edge, it_servers->amount_servers);
                    capacity1.set(edge, 0);
                    capacity2.set(edge, it_servers->amount_servers);
                    cost.set(edge, 0);

                    //Connect this timestep with prior one
                    edge = g.addArc(old_l[k], l_km);
                    capacity.set(edge, it_servers->amount_servers);
                    capacity1.set(edge, it_servers->amount_servers);
                    capacity2.set(edge, it_servers->amount_servers);
                    cost.set(edge, 0);
                }
                SmartDigraph::Arc edge = g.addArc(l_km, l_kma);
                capacity.set(edge, it_servers->amount_servers);
                capacity1.set(edge, it_servers->amount_servers);
                capacity2.set(edge, 0);
                cost.set(edge, it_servers->consumption_rate[k+1][j]);

                edge = g.addArc(sources[j+1], l_kma);
                capacity.set(edge, it_servers->amount_servers);
                capacity1.set(edge, 0);
                capacity2.set(edge, it_servers->amount_servers);
                cost.set(edge, 0);
            }
            old_l = current_l;
        }
        SmartDigraph::Node u_kn = g.addNode();

        for(int i = 0; i != it_servers->transition_costs.size() - 1; ++i) {
            SmartDigraph::Node l_n = g.addNode();
            SmartDigraph::Arc edge = g.addArc(old_l[i], l_n);
            capacity.set(edge, it_servers->amount_servers);
            capacity1.set(edge, it_servers->amount_servers);
            capacity2.set(edge, it_servers->amount_servers);
            cost.set(edge, 0);

            edge = g.addArc(l_n, sinks[sinks.size() - 1]);
            capacity.set(edge, it_servers->amount_servers);
            capacity1.set(edge, 0);
            capacity2.set(edge, it_servers->amount_servers);
            cost.set(edge, 0);
        }

        SmartDigraph::Node l_kn = g.addNode();
        SmartDigraph::Arc edge = g.addArc(u_kn, l_kn);
        capacity.set(edge, it_servers->amount_servers);
        capacity1.set(edge, it_servers->amount_servers);
        capacity2.set(edge, 0);
        cost.set(edge, 0);

        edge = g.addArc(l_kn, sinks[0]);
        capacity.set(edge, it_servers->amount_servers);
        capacity1.set(edge, it_servers->amount_servers);
        capacity2.set(edge, 0);
        cost.set(edge, 0);

        edge = g.addArc(old_u, u_kn);
        capacity.set(edge, it_servers->amount_servers);
        capacity1.set(edge, it_servers->amount_servers);
        capacity2.set(edge, 0);
        cost.set(edge, it_servers->consumption_rate[0][it_servers->consumption_rate.size()-1]);

        edge = g.addArc(old_l[old_l.size()-1], l_kn);
        capacity.set(edge, it_servers->amount_servers);
        capacity1.set(edge, it_servers->amount_servers);
        capacity2.set(edge, it_servers->amount_servers);
        cost.set(edge, 0);

        edge = g.addArc(l_kn, sinks[sinks.size() - 1]);
        capacity.set(edge, it_servers->amount_servers);
        capacity1.set(edge, 0);
        capacity2.set(edge, it_servers->amount_servers);
        cost.set(edge, 0);
    }
}

void print_graph(SmartDigraph &g) {
    ofstream file;
    file.open("dot.gv");
    file << "digraph G {" << std::endl;
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        file  << g.id(g.source(a)) << " -> " << g.id(g.target(a))
              << std::endl;
    }
    file << "}" << std::endl;
}
int main() {
    SmartDigraph g;
    SmartDigraph::NodeMap<long> imbalances1(g);
    SmartDigraph::NodeMap<long> imbalances2(g);
    SmartDigraph::ArcMap<unsigned int> capacity(g);
    SmartDigraph::ArcMap<unsigned int> capacity1(g);
    SmartDigraph::ArcMap<unsigned int> capacity2(g);
    SmartDigraph::ArcMap<unsigned int> cost(g);

    vector<Server> servers;
    vector<unsigned int> demands {1,4};

    vector<vector<unsigned int>> cr {{3,3,3},{2,2,2},{1,1,1}};
    vector<unsigned int> tc {1,2};
    Server server = Server(cr, tc, 3);
    servers.push_back(server);

    vector<vector<unsigned int>> consumption_rate {{5,5,5},{3,3,3},{1,1,1}};
    vector<unsigned int> transition_cost {2,4};
    server = Server(consumption_rate, transition_cost, 3);
    servers.push_back(server);

    generate_graph(g, imbalances1, imbalances2, cost, capacity, capacity1, capacity2, servers, demands);
    print_graph(g);
    return 0;
}