#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <lemon/smart_graph.h>

#include<lemon/glpk.h>
#include<lemon/cplex.h>
#include<lemon/lp.h>

#include<ctime>
#include<sys/timeb.h>
#include<cinttypes>

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
        cost.set(edge, it_servers->consumption_rate[0][it_servers->consumption_rate[0].size()-1]);

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
        file  << g.id(g.source(a)) << " -> " << g.id(g.target(a)) << " [fontcolor=red, label=\"" << g.id(a) << ":" << to_string(cost[a]) << "/" << to_string(capacity[a]) << "/" << to_string(capacity1[a]) << "/" << to_string(capacity2[a]) << "\" ]"
              << std::endl;
    }
    file << "}" << std::endl;
}

void print_graph(const SmartDigraph &g, const SmartDigraph::ArcMap<unsigned long> &capacity, const SmartDigraph::ArcMap<unsigned long> &capacity1,
                 const SmartDigraph::ArcMap<unsigned long> &capacity2, SmartDigraph::ArcMap<unsigned long> &cost, SmartDigraph::NodeMap<long> &imbalances1,
                 SmartDigraph::NodeMap<long> &imbalances2, SmartDigraph::ArcMap<Lp::Col> &f1, SmartDigraph::ArcMap<Lp::Col> &f2, Lp &lp) {
    ofstream file;
    file.open("dot.gv");
    file << "digraph G {" << std::endl;
    for (SmartDigraph::NodeIt a(g); a != INVALID; ++a) {
        file << g.id(a) << " [ xlabel=\"" << to_string(imbalances1[a]) << "/" << to_string(imbalances2[a]) << "\" ]" << std::endl;
    }
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        file  << g.id(g.source(a)) << " -> " << g.id(g.target(a)) << " [fontcolor=red, label=\"" << to_string(cost[a])
              << "|" << to_string(capacity[a]) << "|" << lp.primal(f1[a]) << "/" << to_string(capacity1[a])
              << "|" << lp.primal(f2[a]) << "/" << to_string(capacity2[a]) << "\" ]"
              << std::endl;
    }
    file << "}" << std::endl;
}

void scale_flow(unsigned long amount_servers, const SmartDigraph &g, SmartDigraph::ArcMap<unsigned long> &capacity1,
                SmartDigraph::ArcMap<unsigned long> &cost, SmartDigraph::NodeMap<long> &imbalances1, SmartDigraph::ArcMap<Lp::Col> &f1, Lp &lp) {
    SmartDigraph::ArcMap<double> flow(g);
    for(SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        flow[a] = lp.primal(f1[a]) * amount_servers;
        capacity1[a] *= amount_servers;
    }
    //only nodes with imbalance1 are a_0 and b_0
    //TODO: Potentially erronous
    imbalances1[g.nodeFromId(0)] *= amount_servers;
    imbalances1[g.nodeFromId(1)] *= amount_servers;
}

void round_flow_upper() {

}


double mcmcf(bool is_debug, const SmartDigraph &g, const SmartDigraph::ArcMap<unsigned long> &capacity, const SmartDigraph::ArcMap<unsigned long> &capacity1,
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

        lp.addRow( f1[a] + f2[a] <= capacity[a]);
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
            cerr << "UNBOUNDED" << std::endl;
            break;
        case Lp::ProblemType::OPTIMAL:
            cerr << "OPTIMAL" << std::endl;
                break;
        case Lp::ProblemType::INFEASIBLE:
            cerr << "INFEASIBLE" << std::endl;
            break;
        case Lp::ProblemType::FEASIBLE:
            cerr << "FEASIBLE" << std::endl;
            break;
        case Lp::ProblemType::UNDEFINED:
            cerr << "UNDEFINED" << std::endl;
            break;
    }

    if(is_debug) {
        for(SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
            printf("%d: %f %f\n", g.id(a), lp.primal(f1[a]), lp.primal(f2[a]));
        }
        print_graph(g, capacity, capacity1, capacity2, cost, imbalances1, imbalances2, f1, f2, lp);
    }
    return lp.primal();
}

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

void read_test_file(const string &name, vector<Server> &servers, vector<unsigned long> &demands) {
    string line;
    ifstream file(name);
    if (file.is_open()) {
        getline(file, line);
        vector<string> buffer = split(line, ' ');
        long amount_servers = stol(buffer[0]);
        long amount_demands = stol(buffer[1]);

        getline(file, line);
        buffer = split(line, ' ');
        for (auto it = buffer.begin(); it != buffer.end(); ++it) {
            if (!(*it).empty()) {
                demands.emplace_back(stol(*it));
            }
        }
        for (int i = 0; i != amount_servers; ++i) {
            getline(file, line);
            buffer = split(line, ' ');
            long amount_per_server = stol(buffer[0]);
            long amount_lower = stol(buffer[1]);
            getline(file, line);
            buffer = split(line, ' ');
            vector<unsigned long> tc;
            for(auto it = buffer.begin(); it != buffer.end(); ++it) {
                if(!(*it).empty()) {
                    tc.emplace_back(stol(*it));
                }
            }
            vector<vector<unsigned long>> consumption_rate;
            for (int i = 0; i != amount_lower + 1; ++i) {
                getline(file, line);
                buffer = split(line, ' ');
                vector<unsigned long> cr;
                for (auto it = buffer.begin(); it != buffer.end(); ++it) {
                if (!(*it).empty()) {
                    cr.emplace_back(stol(*it));
                }
            }
                consumption_rate.emplace_back(cr);
            }
            servers.emplace_back(Server(consumption_rate, tc, amount_per_server));
        }
        file.close();
    } else cout << "Unable to open file " << name;
}

uint64_t getTimeNow() {
    timespec ts{};
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t) ts.tv_sec * 1000000LL + (uint64_t) ts.tv_nsec / 1000LL;
}

void benchmark(bool is_debug, int filenumber, uint64_t &generate, uint64_t &flow, string output) {
    srand(time(NULL));
    vector<Server> servers;
    vector<unsigned long> demands;
    read_test_file("tests/test_" + to_string(filenumber), servers, demands);
    ofstream file;
    ofstream data_file;
    file.open(output, std::ios_base::app);
    data_file.open(output + "_data", std::ios_base::app);
    long amount_nodes = 2 * (demands.size() + 1);
    long amount_edges = 0;

    for(auto it_servers = servers.begin(); it_servers != servers.end(); ++it_servers) {
        amount_nodes += (1 + 2 * it_servers->transition_costs.size()) * demands.size() + 1 + it_servers->transition_costs.size();
        amount_edges += 5 + 4 * it_servers->transition_costs.size() + (1 + 6 * it_servers->transition_costs.size())
 * (demands.size() - 1);
    }

    SmartDigraph g;
    SmartDigraph::NodeMap<long> imbalances1(g);
    SmartDigraph::NodeMap<long> imbalances2(g);
    SmartDigraph::ArcMap<unsigned long> capacity(g);
    SmartDigraph::ArcMap<unsigned long> capacity1(g);
    SmartDigraph::ArcMap<unsigned long> capacity2(g);
    SmartDigraph::ArcMap<unsigned long> cost(g);

    uint64_t start = getTimeNow();
    generate_graph(g, imbalances1, imbalances2, cost, capacity, capacity1, capacity2, servers, demands);
    uint64_t start_flow = getTimeNow();

    double min_cost = mcmcf(is_debug, g, capacity, capacity1, capacity2, cost, imbalances1, imbalances2);
    uint64_t end = getTimeNow();

    cout << "Minimal Cost: " << min_cost << std::endl;

    uint64_t generate_bench = (start_flow - start);
    uint64_t flow_bench = (end - start_flow);
    uint64_t bench = generate_bench + flow_bench;
    generate += generate_bench;
    flow += flow_bench;
    printf("Generating the graph took %" PRIu64 " nanoseconds.\n", generate_bench);
    printf("Executing mcmcf took %" PRIu64 " nanoseconds.\n", flow_bench);
    printf("Overall time required was %" PRIu64 " nanoseconds.\n", bench);
    file << (filenumber + 1) << " & " << amount_nodes << " & " << amount_edges << " & ";
    file << "\\SI{" << generate_bench << "}{\\nano\\second} & " << "\\SI{" << flow_bench << "}{\\nano\\second} & "
         << "\\SI{" << bench << "}{\\nano\\second}\\\\" << std::endl;
    file << "\\hline" << std::endl;
    data_file << generate_bench << " & " << flow_bench << " & " << bench << "\\\\" << std::endl;
    file.close();
}

int main(int argc, char * argv[]) {
    /*
    string testnumber = argv[1];
    SmartDigraph g;
    SmartDigraph::NodeMap<long> imbalances1(g);
    SmartDigraph::NodeMap<long> imbalances2(g);
    SmartDigraph::ArcMap<unsigned long> capacity(g);
    SmartDigraph::ArcMap<unsigned long> capacity1(g);
    SmartDigraph::ArcMap<unsigned long> capacity2(g);
    SmartDigraph::ArcMap<unsigned long> cost(g);

    vector<Server> servers;
    vector<unsigned long> demands {1,4};

    vector<vector<unsigned long>> cr {{3,3},{2,2},{1,1}};
    vector<unsigned long> tc {1,2};
    Server server = Server(cr, tc, 3);
    servers.push_back(server);

    vector<vector<unsigned long>> consumption_rate {{5,5},{3,3},{1,1}};
    vector<unsigned long> transition_cost {2,4};
    server = Server(consumption_rate, transition_cost, 3);
    servers.push_back(server);
    read_test_file("tests/test_" + testnumber, servers, demands);
    generate_graph(g, imbalances1, imbalances2, cost, capacity, capacity1, capacity2, servers, demands);
    print_graph(g, capacity, capacity1, capacity2, cost, imbalances1, imbalances2);
    double min_cost = mcmcf(g, capacity, capacity1, capacity2, cost, imbalances1, imbalances2);
    cout << "Minimal Cost: " << min_cost << std::endl;
    */
    const int amount_tests = stoi(argv[1]);

    uint64_t generate = 0;
    uint64_t flow = 0;
    for(int i = 0; i != amount_tests; ++i) {
        printf("running test %d\n", i);
        benchmark(true, i, generate, flow, "result");
    }

    printf("Ran %d tests.\n", amount_tests);
    printf("Overall time required for generating the graph: %" PRIu64 " nanoseconds\n", generate);
    printf("Overall time required for executing the simplex algorithm: %" PRIu64 " nanoseconds\n", flow);
    printf("Overall time required: %" PRIu64 " nanoseconds\n", generate + flow);
    return 0;
}