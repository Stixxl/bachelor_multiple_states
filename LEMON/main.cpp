#include <iostream>
#include <fstream>
#include <vector>
#include <lemon/smart_graph.h>

#include<lemon/glpk.h>
#include<lemon/lp.h>

using namespace std;
using namespace lemon;

struct Server {
    vector<vector<unsigned long>> consumption_rate;
    vector<unsigned long> transition_costs;
    unsigned long amount_servers;

    Server(vector<vector<unsigned long>> cr, vector<unsigned long> tc, unsigned long amount) {
        consumption_rate = move(cr);
        transition_costs = move(tc);
        amount_servers = amount;
    }
};

void generate_graph(SmartDigraph &g, SmartDigraph::NodeMap<long> &imbalances1, SmartDigraph::NodeMap<long> &imbalances2,
                    SmartDigraph::ArcMap<unsigned long> &cost, SmartDigraph::ArcMap<unsigned long> &capacity,
                    SmartDigraph::ArcMap<unsigned long> &capacity1, SmartDigraph::ArcMap<unsigned long> &capacity2,
                    vector<Server> &servers, vector<unsigned long> &demands) {
    unsigned long amount_servers = 0;
    unsigned long d_k = 0;
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
        for(long j = 0; j != demands.size(); ++j) {
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
            for(long k = 0; k != it_servers->transition_costs.size(); ++k) {
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

        for(long i = 0; i != it_servers->transition_costs.size() - 1; ++i) {
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

void print_graph(const SmartDigraph &g, const SmartDigraph::ArcMap<unsigned long> &capacity, const SmartDigraph::ArcMap<unsigned long> &capacity1,
                 const SmartDigraph::ArcMap<unsigned long> &capacity2, SmartDigraph::ArcMap<unsigned long> &cost, SmartDigraph::NodeMap<long> &imbalances1,
                 SmartDigraph::NodeMap<long> &imbalances2) {
    ofstream file;
    file.open("dot.gv");
    file << "digraph G {" << std::endl;
    for (SmartDigraph::NodeIt a(g); a != INVALID; ++a) {
        file << g.id(a) << " [ xlabel=\"" << to_string(imbalances1[a]) << "/" << to_string(imbalances2[a]) << "\" ]" << std::endl;
    }
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        file  << g.id(g.source(a)) << " -> " << g.id(g.target(a)) << " [fontcolor=red, label=\"" << to_string(cost[a]) << "/" << to_string(capacity[a]) << "/" << to_string(capacity1[a]) << "/" << to_string(capacity2[a]) << "\" ]"
              << std::endl;
    }
    file << "}" << std::endl;
}


double mcmcf(const SmartDigraph &g, const SmartDigraph::ArcMap<unsigned long> &capacity, const SmartDigraph::ArcMap<unsigned long> &capacity1,
const SmartDigraph::ArcMap<unsigned long> &capacity2, SmartDigraph::ArcMap<unsigned long> &cost, SmartDigraph::NodeMap<long> &imbalances1,
             SmartDigraph::NodeMap<long> &imbalances2)
{
    // Create an instance of the default LP solver
    // Add a column to the problem for each arc
    Lp lp;
    SmartDigraph::ArcMap<Lp::Col> f1(g);
    SmartDigraph::ArcMap<Lp::Col> f2(g);
    lp.addColSet(f1);
    lp.addColSet(f2);
    // Capacity constraints
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        lp.colLowerBound(f1[a], 0);
        lp.colUpperBound(f1[a], capacity1[a]);

        lp.colLowerBound(f2[a], 0);
        lp.colUpperBound(f2[a], capacity2[a]);

        lp.addRow(f1[a] + f2[a] <= capacity[a]);
        //lp.addRow(0 <= f1[a] + f2[a]); implicitly set since both f1[a] and f2[a] are required to be >= 0
    }
    // Flow conservation constraints
    for (SmartDigraph::NodeIt n(g); n != INVALID; ++n) {
        Lp::Expr e;
        Lp::Expr e1;
        Lp::Expr e2;
        for (SmartDigraph::OutArcIt a(g, n); a != INVALID; ++a)
        {
            e1 += f1[a];
            e2 += f2[a];
            e += f1[a] + f2[a];
        }
        for (SmartDigraph::InArcIt a(g, n); a != INVALID; ++a) {
            e1 -= f1[a];
            e2 -= f2[a];
            e -= f1[a] + f2[a];
        }
        lp.addRow(e == imbalances1[n] + imbalances2[n]);
        lp.addRow(e1 == imbalances1[n]);
        lp.addRow(e2 == imbalances2[n]);
    }
    // Objective function
    Lp::Expr o;
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        o += (f1[a] + f2[a]) * cost[a];
    }
    lp.min();
    lp.obj(o);
    // Solve the LP problem
    lp.solve();
    switch(lp.primalType()) {
        case Lp::ProblemType::UNBOUNDED:
            cout << "UNBOUNDED" << std::endl;
            break;
        case Lp::ProblemType::OPTIMAL:
            cout << "OPTIMAL" << std::endl;
                break;
        case Lp::ProblemType::INFEASIBLE:
            cout << "INFEASIBLE" << std::endl;
            break;
        case Lp::ProblemType::FEASIBLE:
            cout << "FEASIBLE" << std::endl;
            break;
        case Lp::ProblemType::UNDEFINED:
            cout << "UNDEFINED" << std::endl;
            break;
    }

    return lp.primal();
}
int main() {
    SmartDigraph g;
    SmartDigraph::NodeMap<long> imbalances1(g);
    SmartDigraph::NodeMap<long> imbalances2(g);
    SmartDigraph::ArcMap<unsigned long> capacity(g);
    SmartDigraph::ArcMap<unsigned long> capacity1(g);
    SmartDigraph::ArcMap<unsigned long> capacity2(g);
    SmartDigraph::ArcMap<unsigned long> cost(g);

    vector<Server> servers;
    vector<unsigned long> demands {1,4};

    vector<vector<unsigned long>> cr {{3,3,3},{2,2,2},{1,1,1}};
    vector<unsigned long> tc {1,2};
    Server server = Server(cr, tc, 3);
    servers.push_back(server);

    vector<vector<unsigned long>> consumption_rate {{5,5,5},{3,3,3},{1,1,1}};
    vector<unsigned long> transition_cost {2,4};
    server = Server(consumption_rate, transition_cost, 3);
    servers.push_back(server);

    generate_graph(g, imbalances1, imbalances2, cost, capacity, capacity1, capacity2, servers, demands);
    print_graph(g, capacity, capacity1, capacity2, cost, imbalances1, imbalances2);
    double min_cost = mcmcf(g, capacity, capacity1, capacity2, cost, imbalances1, imbalances2);
    cout << "Minimal Cost: " << min_cost << std::endl;
    return 0;
}