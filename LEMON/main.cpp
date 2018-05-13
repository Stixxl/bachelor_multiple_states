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
#include <chrono>
#include<cinttypes>
#include<limits>

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCDFAInspection"
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
    for (auto it_server = servers.begin(); it_server != servers.end(); ++it_server) {
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

    for (auto it_demands = demands.begin(); it_demands != demands.end(); ++it_demands) {
        SmartDigraph::Node a_k = g.addNode();
        imbalances2.set(a_k, d_k + *it_demands);
        sources.push_back(a_k);

        SmartDigraph::Node b_k = g.addNode();
        imbalances2.set(b_k, -1 * (d_k + *it_demands));
        sinks.push_back(b_k);
    }

    SmartDigraph::Node old_u;
    vector<SmartDigraph::Node> old_l;
    for (auto it_servers = servers.begin(); it_servers != servers.end(); ++it_servers) {
        for (long j = 0; j != demands.size(); ++j) {
            SmartDigraph::Node u_k = g.addNode();
            vector<SmartDigraph::Node> current_l;
            if (j != 0) {
                //Connect this timestep with prior one
                SmartDigraph::Arc edge = g.addArc(old_u, u_k);
                capacity.set(edge, it_servers->amount_servers);
                capacity1.set(edge, it_servers->amount_servers);
                capacity2.set(edge, 0);
                cost.set(edge, it_servers->consumption_rate[0][j - 1]);
            }
            old_u = u_k;
            for (long k = 0; k != it_servers->transition_costs.size(); ++k) {
                SmartDigraph::Node l_km = g.addNode();
                SmartDigraph::Node l_kma = g.addNode();
                current_l.push_back(l_kma);
                if (j == 0 && k == it_servers->transition_costs.size() - 1) {
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
                if (j != 0) {
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
                cost.set(edge, it_servers->consumption_rate[k + 1][j]);

                edge = g.addArc(sources[j + 1], l_kma);
                capacity.set(edge, it_servers->amount_servers);
                capacity1.set(edge, 0);
                capacity2.set(edge, it_servers->amount_servers);
                cost.set(edge, 0);
            }
            old_l = current_l;
        }
        SmartDigraph::Node u_kn = g.addNode();

        for (long i = 0; i != it_servers->transition_costs.size() - 1; ++i) {
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
        cost.set(edge, it_servers->consumption_rate[0][it_servers->consumption_rate[0].size() - 1]);

        edge = g.addArc(old_l[old_l.size() - 1], l_kn);
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

void print_graph(const SmartDigraph &g, const SmartDigraph::ArcMap<unsigned long> &capacity,
                 const SmartDigraph::ArcMap<unsigned long> &capacity1,
                 const SmartDigraph::ArcMap<unsigned long> &capacity2, SmartDigraph::ArcMap<unsigned long> &cost,
                 SmartDigraph::NodeMap<long> &imbalances1,
                 SmartDigraph::NodeMap<long> &imbalances2) {
    ofstream file;
    file.open("results/dot.gv");
    file << "digraph G {" << std::endl;
    for (SmartDigraph::NodeIt a(g); a != INVALID; ++a) {
        file << g.id(a) << " [ xlabel=\"" << to_string(imbalances1[a]) << "/" << to_string(imbalances2[a]) << "\" ]"
             << std::endl;
    }
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        file << g.id(g.source(a)) << " -> " << g.id(g.target(a)) << " [fontcolor=red, label=\"" << g.id(a) << ":"
             << to_string(cost[a]) << "/" << to_string(capacity[a]) << "/" << to_string(capacity1[a]) << "/"
             << to_string(capacity2[a]) << "\" ]"
             << std::endl;
    }
    file << "}" << std::endl;
}

void print_graph(const SmartDigraph &g, const SmartDigraph::ArcMap<unsigned long> &capacity,
                 const SmartDigraph::ArcMap<unsigned long> &capacity1,
                 const SmartDigraph::ArcMap<unsigned long> &capacity2, SmartDigraph::ArcMap<unsigned long> &cost,
                 SmartDigraph::NodeMap<long> &imbalances1,
                 SmartDigraph::NodeMap<long> &imbalances2, SmartDigraph::ArcMap<Lp::Col> &f1,
                 SmartDigraph::ArcMap<Lp::Col> &f2, Lp &lp) {
    ofstream file;
    file.open("dot.gv");
    file << "digraph G {" << std::endl;
    for (SmartDigraph::NodeIt a(g); a != INVALID; ++a) {
        file << g.id(a) << " [ xlabel=\"" << to_string(imbalances1[a]) << "/" << to_string(imbalances2[a]) << "\" ]"
             << std::endl;
    }
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        file << g.id(g.source(a)) << " -> " << g.id(g.target(a)) << " [fontcolor=red, label=\"" << to_string(cost[a])
             << "|" << to_string(capacity[a]) << "|" << lp.primal(f1[a]) << "/" << to_string(capacity1[a])
             << "|" << lp.primal(f2[a]) << "/" << to_string(capacity2[a]) << "\" ]"
             << std::endl;
    }
    file << "}" << std::endl;
}

void scale_flow(unsigned long amount_servers, const SmartDigraph &g, SmartDigraph::ArcMap<unsigned long> &capacity1,
                SmartDigraph::ArcMap<double> &flow, SmartDigraph::NodeMap<long> &imbalances1,
                SmartDigraph::ArcMap<Lp::Col> &f1, Lp &lp) {
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        flow[a] = lp.primal(f1[a]) * amount_servers;
        capacity1[a] *= amount_servers;
    }
    //only nodes with imbalance1 are a_0 and b_0
    imbalances1[g.nodeFromId(0)] *= amount_servers;
    imbalances1[g.nodeFromId(1)] *= amount_servers;
}

void round_flow_valley(ArcLookUp<SmartDigraph> &ae, const vector<Server> &servers, const unsigned long amnt_timesteps,
                       const SmartDigraph &g,
                       SmartDigraph::ArcMap<double> &flow) {
    //last step does not have to be considered since the flow there is always decreasing
    unsigned long vertices_counter = 2 * (amnt_timesteps + 1);
    unsigned long server_counter = 0;
    double prior_flow = 0.0;
    bool is_decreasing = false;
    unsigned long valley_start = 0;
    unsigned long path_length = 0;
    double overflow = 0.0;
    while (server_counter < servers.size()) {
        Server current = servers[server_counter];
        unsigned long sigma = current.transition_costs.size();
        unsigned long u_k = vertices_counter;
        unsigned long u_k1 = vertices_counter + 1 + 2 * sigma;
        SmartDigraph::Arc edge = ae(g.nodeFromId(u_k), g.nodeFromId(u_k1));
        if (edge == INVALID) {
            vertices_counter += sigma + 1;
            server_counter++;
            is_decreasing = false;
            continue;
        } else {
            //found end of valley
            if (flow[edge] > prior_flow && is_decreasing) {
                //find smallest lower path that routes flow
                unsigned long l_k = valley_start + 1;
                for (int i = 0; i != sigma; ++i) {
                    SmartDigraph::Arc edge = ae(g.nodeFromId(u_k), g.nodeFromId(l_k));
                    if (flow[edge] > 0.0) {
                        break;
                    } else {
                        l_k += 2;
                    }
                }
                //route flow to lower path
                SmartDigraph::Arc power_down = ae(g.nodeFromId(valley_start), g.nodeFromId(l_k));
                flow[power_down] += overflow;
                u_k = valley_start;
                u_k1 = valley_start + 1 + 2 * sigma;
                unsigned long l_ak = l_k + 1;
                unsigned long l_k1 = l_k + 2 * sigma + 1;
                for (int i = 0; i != path_length; ++i) {
                    //remove flow from upper path
                    SmartDigraph::Arc upper = ae(g.nodeFromId(u_k), g.nodeFromId(u_k1));
                    flow[upper] -= overflow;
                    //route flow over lower paths
                    SmartDigraph::Arc lower = ae(g.nodeFromId(l_k), g.nodeFromId(l_ak));
                    SmartDigraph::Arc lower1 = ae(g.nodeFromId(l_ak), g.nodeFromId(l_k1));
                    flow[lower] += overflow;
                    flow[lower1] += overflow;
                    l_k = l_k1;
                    l_ak = l_k + 1;
                    l_k1 += 1 + 2 * sigma;
                    u_k = u_k1;
                    u_k1 += 1 + 2 * sigma;
                }
                //route flow back to upper path
                SmartDigraph::Arc power_up = ae(g.nodeFromId(l_k), g.nodeFromId(u_k));
                flow[power_up] += overflow;
                is_decreasing = false;
                vertices_counter = u_k1;
                edge = ae(g.nodeFromId(u_k), g.nodeFromId(u_k1));
                //found potential start of valley
            } else {
                vertices_counter += 1 + 2 * sigma;
                path_length++;
            }
        }
        prior_flow = flow[edge];
    }
}

void round_flow_inc(ArcLookUp<SmartDigraph> &ae, const vector<Server> &servers, const unsigned long amnt_timesteps,
                    const SmartDigraph &g,
                    SmartDigraph::ArcMap<double> &flow) {
    //last step does not have to be considered since the flow there is always decreasing
    unsigned long vertices_counter = 2 * (amnt_timesteps + 1);
    unsigned long server_counter = 0;
    double prior_flow = 0.0;
    double overflow;
    bool is_decreasing = true;
    while (server_counter < servers.size()) {
        unsigned long sigma = servers[server_counter].transition_costs.size();
        unsigned long u_k = vertices_counter;
        unsigned long u_k1 = vertices_counter + 1 + 2 * sigma;
        SmartDigraph::Arc edge = ae(g.nodeFromId(u_k), g.nodeFromId(u_k1));
        if (edge == INVALID) {
            vertices_counter += sigma + 1;
            is_decreasing = true;
            server_counter++;
            continue;
        } else {
            if (flow[edge] > prior_flow & is_decreasing) {
                //FIXME Check if it works
                //following 3 comments also count for decreasing and peaks
                //FIXME u_1 also counts as end of valley
                //FIXME vertex at which flow decreases; prior edge will not be updated
                //FIXME per update procedure: ship overflow until flow is integral

                is_decreasing = false;
                //FIXME Check whether algorithm updates correct vertices in between
                while (!is_decreasing) {
                    u_k = vertices_counter;
                    u_k1 = vertices_counter + 1 + 2 * sigma;
                    SmartDigraph::Arc upper = ae(g.nodeFromId(u_k), g.nodeFromId(u_k1));
                    prior_flow = flow[upper];
                    //exit if (u_{k''-1}, u_{k''}) is reached

                    double delta = flow[upper] - floor(flow[upper]);
                    while (delta > 0.0) {
                        //find largest lower path that routes flow
                        unsigned long l_k = u_k + 2 * sigma - 1;
                        for (int i = 0; i != sigma; ++i) {
                            SmartDigraph::Arc edge = ae(g.nodeFromId(l_k), g.nodeFromId(u_k));
                            if (flow[edge] > 0.0) {
                                break;
                            } else {
                                l_k -= 2;
                            }
                        }
                        unsigned long l_ak = l_k + 1;
                        unsigned long l_k1 = l_k + 2 * sigma + 1;
                        SmartDigraph::Arc power_up = ae(g.nodeFromId(l_k), g.nodeFromId(u_k));
                        overflow = min(delta, flow[power_up]);
                        flow[upper] -= overflow;
                        //route flow to lower path
                        flow[power_up] -= overflow;
                        //route flow over lower paths
                        SmartDigraph::Arc lower = ae(g.nodeFromId(l_k), g.nodeFromId(l_ak));
                        SmartDigraph::Arc lower1 = ae(g.nodeFromId(l_ak), g.nodeFromId(l_k1));
                        SmartDigraph::Arc power_up1 = ae(g.nodeFromId(l_k1), g.nodeFromId(u_k1));
                        flow[lower] += overflow;
                        flow[lower1] += overflow;
                        flow[power_up1] += overflow;
                        delta -= overflow;
                    }
                    vertices_counter += 1 + 2 * sigma;
                    SmartDigraph::Arc next_upper = ae(g.nodeFromId(vertices_counter),
                                                      g.nodeFromId(vertices_counter + 1 + 2 * sigma));
                    if (next_upper == INVALID) {
                        if (flow[upper] > 0) {
                            is_decreasing = true;
                        }
                    } else {
                        is_decreasing = flow[next_upper] < flow[upper];
                    }
                }
            } else if (flow[edge] < prior_flow) {
                is_decreasing = true;
                prior_flow = flow[edge];
                vertices_counter += 1 + 2 * sigma;//set counter to next u_{i,k +1 }
            } else {
                prior_flow = flow[edge];
                vertices_counter += 1 + 2 * sigma;//set counter to next u_{i,k +1 }
            }
        }
    }
}

void round_flow_dec(ArcLookUp<SmartDigraph> &ae, const vector<Server> &servers, const unsigned long amnt_timesteps,
                    const SmartDigraph &g,
                    SmartDigraph::ArcMap<double> &flow) {
    //last step does not have to be considered since the flow there is always decreasing
    unsigned long vertices_counter = 2 * (amnt_timesteps + 1) + 1 + 2 * servers[0].transition_costs.size();
    unsigned long u_n = vertices_counter;
    unsigned long server_counter = 0;
    bool is_increasing = true;
    double overflow;
    while (server_counter < servers.size()) {
        unsigned long sigma = servers[server_counter].transition_costs.size();
        unsigned long u_k_pred = vertices_counter - (1 + 2 * sigma);
        unsigned long u_k = vertices_counter;
        unsigned long u_k1 = vertices_counter + 1 + 2 * sigma;
        SmartDigraph::Arc edge = ae(g.nodeFromId(u_k_pred), g.nodeFromId(u_k));
        SmartDigraph::Arc edge_next = ae(g.nodeFromId(u_k), g.nodeFromId(u_k1));
        if (edge == INVALID) {
            vertices_counter = u_n + sigma + 1 + 1 + 2 * sigma;
            u_n = vertices_counter;
            server_counter++;
            is_increasing = true;
            continue;
        } else {
            if (flow[edge] > flow[edge_next] && is_increasing) {
                while (!is_increasing) {
                    u_k = vertices_counter;
                    u_k_pred = vertices_counter - 1 - 2 * sigma;
                    u_k1 = vertices_counter + 1 + 2 * sigma;
                    SmartDigraph::Arc prior_upper = ae(g.nodeFromId(u_k_pred), g.nodeFromId(u_k));
                    SmartDigraph::Arc upper = ae(g.nodeFromId(u_k), g.nodeFromId(u_k1));
                    double delta = flow[upper] - floor(flow[upper]);
                    while (delta > 0.0) {
                        unsigned long l_k = u_k + 2 * (sigma - 1);
                        for (int i = 0; i != sigma; ++i) {
                            SmartDigraph::Arc edge = ae(g.nodeFromId(u_k), g.nodeFromId(l_k));
                            if (flow[edge] > 0.0) {
                                break;
                            } else {
                                l_k -= 2;
                            }
                        }
                        unsigned long l_k_pred = l_k - 2 * sigma;
                        unsigned long l_ak_pred = l_k_pred + 1;
                        //remove flow from upper path
                        SmartDigraph::Arc power_down = ae(g.nodeFromId(u_k), g.nodeFromId(l_k));
                        overflow = min(delta, flow[power_down]);
                        flow[upper] -= overflow;
                        //route flow to lower path
                        flow[power_down] -= overflow;
                        //route flow over lower paths
                        power_down = ae(g.nodeFromId(u_k_pred), g.nodeFromId(l_k_pred));
                        flow[power_down] += overflow;
                        SmartDigraph::Arc lower = ae(g.nodeFromId(l_k_pred), g.nodeFromId(l_ak_pred));
                        SmartDigraph::Arc lower1 = ae(g.nodeFromId(l_ak_pred), g.nodeFromId(l_k));
                        flow[lower] += overflow;
                        flow[lower1] += overflow;

                        delta -= overflow;
                    }
                    if (prior_upper == INVALID) {
                        is_increasing = flow[upper] > 0.0;
                    } else {
                        is_increasing = flow[prior_upper] < flow[upper];
                    }
                }
                vertices_counter = u_k_pred;

            } else if (flow[edge] < flow[edge_next]) {
                vertices_counter -= 1 + 2 * sigma;//set counter to next u_{i,k +1 }
                is_increasing = true;
            } else {
                vertices_counter -= 1 + 2 * sigma;//set counter to next u_{i,k +1 }
            }
        }
    }
}

void round_peaks(ArcLookUp<SmartDigraph> &ae, const vector<Server> &servers, const unsigned long amnt_timesteps,
                 const SmartDigraph &g,
                 SmartDigraph::ArcMap<double> &flow) {
    unsigned long vertices_counter = amnt_timesteps * 2;
    unsigned long server_counter = 0;
    double prior_flow = 0.0;
    bool is_increasing = false;

    while (server_counter < servers.size()) {
        unsigned long sigma = servers[server_counter].transition_costs.size();
        unsigned long u_k = vertices_counter;
        unsigned long u_k1 = vertices_counter + 1 + 2 * sigma;

        SmartDigraph::Arc upper = ae(g.nodeFromId(u_k), g.nodeFromId(u_k1));
        if (upper == INVALID) {
            vertices_counter += 1 + sigma;
            server_counter++;
            is_increasing = false;
            prior_flow = 0.0;
            continue;
        }
        if (prior_flow < flow[upper]) {
            is_increasing = true;
        } else if ((prior_flow > flow[upper] ||
                    (ae(g.nodeFromId(u_k1), g.nodeFromId(u_k1 + 1 + 2 * sigma)) == INVALID && flow[upper] > 0)) &&
                   is_increasing && flow[upper] - floor(flow[upper] > 0.0)) {
            unsigned long l_k = u_k + 2 * sigma - 1;
            int j1 = sigma;
            for (int i = 0; i != sigma; ++i) {
                SmartDigraph::Arc edge = ae(g.nodeFromId(l_k), g.nodeFromId(u_k));
                if (flow[edge] > 0.0) {
                    break;
                } else {
                    l_k -= 2;
                    j1--;
                }
            }
            unsigned long l_k1 = u_k1 + 2 * sigma +1;
            int j2 = sigma;
            for (int i = 0; i != sigma; ++i) {
                SmartDigraph::Arc edge = ae(g.nodeFromId(u_k), g.nodeFromId(l_k));
                if (flow[edge] > 0.0) {
                    break;
                } else {
                    l_k1 -= 2;
                    j2--;
                }
            }
            SmartDigraph::Arc power_up = ae(g.nodeFromId(l_k), g.nodeFromId(u_k));
            SmartDigraph::Arc power_down = ae(g.nodeFromId(u_k1), g.nodeFromId(l_k1));
            if (j1 == j2) {
                double overflow = min(flow[upper] - floor(flow[upper]), flow[power_up]);
                overflow = min(overflow, flow[power_down]);
                flow[upper] -= overflow;
                flow[power_up] -= overflow;
                flow[power_down] -= overflow;

                SmartDigraph::Arc lower = ae(g.nodeFromId(l_k), g.nodeFromId(l_k + 1));
                SmartDigraph::Arc lower1 = ae(g.nodeFromId(l_k), g.nodeFromId(l_k1));
                flow[lower] += overflow;
                flow[lower1] += overflow;
            } else {
                double overflow = flow[upper] - floor(flow[upper]);
                unsigned long l_k_pred = l_k - 2 * sigma - 1;
                SmartDigraph::Arc edge = ae(g.nodeFromId(l_k_pred), g.nodeFromId(l_k));
                while (edge != INVALID && flow[edge] > 0.0) {
                    overflow = min(overflow, flow[edge]);
                    l_k = l_k_pred;
                    l_k_pred -= 2 * sigma - 1;
                    edge = ae(g.nodeFromId(l_k_pred), g.nodeFromId(l_k));
                }
                unsigned long l_k1_old = l_k1;
                unsigned long l_k2 = l_k1 + 2 * sigma + 1;
                edge = ae(g.nodeFromId(l_k1), g.nodeFromId(l_k2));
                while (edge != INVALID && flow[edge] > 0.0) {
                    overflow = min(overflow, flow[edge]);
                    l_k1 = l_k2;
                    l_k2 += 2 * sigma + 1;
                    edge = ae(g.nodeFromId(l_k1), g.nodeFromId(l_k2));
                }
                flow[upper] -= overflow;
                if (j1 < j2) {
                    unsigned long start = l_k;
                    unsigned long start_lower = l_k + 2 * (j2 - j1);
                    unsigned long u_start = l_k - 2 * j1 - 1;
                    SmartDigraph::Arc edge = ae(g.nodeFromId(u_start), g.nodeFromId(start));
                    flow[edge] -= overflow;
                    edge = ae(g.nodeFromId(u_start), g.nodeFromId(start_lower));
                    flow[edge] += overflow;
                    while (start_lower != l_k1_old) {
                        unsigned long start_next = start + 2 * sigma + 1;
                        unsigned long start_lower_next = start_lower + 2 * sigma + 1;
                        SmartDigraph::Arc edge = ae(g.nodeFromId(start), g.nodeFromId(start_next));
                        flow[edge] -= overflow;
                        edge = ae(g.nodeFromId(start_lower), g.nodeFromId(start_lower_next));
                        flow[edge] += overflow;
                        start = start_next;
                        start_lower = start_lower_next;
                    }
                    unsigned long u_k_end = start - j1 * 2 - 1;
                    edge = ae(g.nodeFromId(start), g.nodeFromId(u_k_end));
                    flow[edge] -= overflow;
                    edge = ae(g.nodeFromId(start_lower), g.nodeFromId(start_lower + 2 * sigma + 1));
                    flow[edge] += overflow;
                    edge = ae(g.nodeFromId(u_k1), g.nodeFromId(start_lower + 2 * sigma + 1));
                    flow[edge] -= overflow;
                } else {
                    unsigned long start = u_k1 + j2 * 2 - 1;
                    unsigned long start_lower = u_k1 + j1 * 2 - 1;
                    SmartDigraph::Arc edge = ae(g.nodeFromId(u_k1), g.nodeFromId(start_lower));
                    flow[edge] -= overflow;
                    edge = ae(g.nodeFromId(start_lower), g.nodeFromId(start_lower + 2 * sigma + 1));
                    flow[edge] += overflow;
                    start += 2 * sigma + 1;
                    start_lower += 2 * sigma + 1;
                    while (start != l_k1) {
                        unsigned long start_next = start + 2 * sigma + 1;
                        unsigned long start_lower_next = start_lower + 2 * sigma + 1;
                        SmartDigraph::Arc edge = ae(g.nodeFromId(start), g.nodeFromId(start_next));
                        flow[edge] -= overflow;
                        edge = ae(g.nodeFromId(start_lower), g.nodeFromId(start_lower_next));
                        flow[edge] += overflow;
                        start = start_next;
                        start_lower = start_lower_next;
                    }
                    unsigned long u_k_end = start - j2 * 2 - 1;
                    edge = ae(g.nodeFromId(start), g.nodeFromId(u_k_end));
                    flow[edge] -= overflow;
                    edge = ae(g.nodeFromId(start_lower), g.nodeFromId(u_k_end));
                    flow[edge] += overflow;
                }
            }
        } else {
            is_increasing = false;
        }
        prior_flow = flow[upper];
        vertices_counter += 1 + 2 * sigma;
    }
}

void reduce_to_m(ArcLookUp<SmartDigraph> &ae, const vector<Server> &servers, const unsigned long amnt_timesteps,
                 const SmartDigraph &g,
                 SmartDigraph::ArcMap<unsigned long> &capacity1, SmartDigraph::ArcMap<double> &flow) {
    //FIXME most certainly requires fixing
    unsigned long vertices_counter = 2 * amnt_timesteps;
    unsigned long server_counter = 0;

    while (server_counter < servers.size()) {
        unsigned long sigma = servers[server_counter].transition_costs.size();
        unsigned long u_k = vertices_counter;
        unsigned long u_k1 = vertices_counter + 1 + 2 * sigma;
        SmartDigraph::Arc upper = ae(g.nodeFromId(u_k), g.nodeFromId(u_k1));

        if (upper == INVALID) {
            server_counter++;
            vertices_counter += 1 + sigma;
        } else {
            if (flow[upper] > (double) servers[server_counter].amount_servers) {
                unsigned long u_k2 = u_k1 + 1 + 2 * sigma;
                SmartDigraph::Arc next_upper = ae(g.nodeFromId(u_k1), g.nodeFromId(u_k2));
                if (flow[upper] > flow[next_upper] ||
                    (ae(g.nodeFromId(u_k2), g.nodeFromId(u_k2 + 1 + 2 * sigma)) == INVALID &&
                     flow[next_upper] > 0)) {
                    cout << "handling peak" << std::endl;
                    unsigned long l_k = u_k + 2 * sigma - 1;
                    int j1 = sigma-1;
                    for (int i = 0; i != sigma - 2; ++i) {
                        SmartDigraph::Arc edge = ae(g.nodeFromId(l_k), g.nodeFromId(u_k));
                        if (flow[edge] > 0.0) {
                            break;
                        } else {
                            l_k -= 2;
                            j1--;
                        }
                    }
                    unsigned long l_k1 = u_k1 + 2 * sigma - 1;
                    int j2 = sigma - 1;
                    for (int i = 0; i != sigma - 2; ++i) {
                        SmartDigraph::Arc edge = ae(g.nodeFromId(u_k), g.nodeFromId(l_k));
                        if (flow[edge] > 0.0) {
                            break;
                        } else {
                            l_k1 -= 2;
                            j2--;
                        }
                    }
                    double delta;
                    if (next_upper == INVALID) {
                        delta = flow[upper] - (double) servers[server_counter].amount_servers;
                    } else {
                        delta = flow[upper] -
                                max(flow[next_upper], (double) servers[server_counter].amount_servers);
                    }
                    while (flow[upper] > (double) servers[server_counter].amount_servers) {

                        SmartDigraph::Arc power_up = ae(g.nodeFromId(l_k), g.nodeFromId(u_k));
                        SmartDigraph::Arc power_down = ae(g.nodeFromId(u_k1), g.nodeFromId(l_k1));
                        if (j1 == j2) {
                            cout << "j1: " << j1 << ", l_k: " << l_k << ", u_k: " << u_k << std::endl;
                            double overflow = min(delta, flow[power_up]);
                            overflow = min(overflow, flow[power_down]);
                            flow[upper] -= overflow;
                            flow[power_up] -= overflow;
                            flow[power_down] -= overflow;

                            SmartDigraph::Arc lower = ae(g.nodeFromId(l_k), g.nodeFromId(l_k + 1));
                            SmartDigraph::Arc lower1 = ae(g.nodeFromId(l_k), g.nodeFromId(l_k1));
                            flow[lower] += overflow;
                            flow[lower1] += overflow;
                        } else {
                            double overflow = delta;
                            unsigned long l_k_pred = l_k - 2 * sigma - 1;
                            SmartDigraph::Arc edge = ae(g.nodeFromId(l_k_pred), g.nodeFromId(l_k));
                            while (edge != INVALID && flow[edge] > 0.0) {
                                overflow = min(overflow, flow[edge]);
                                l_k = l_k_pred;
                                l_k_pred -= 2 * sigma - 1;
                                edge = ae(g.nodeFromId(l_k_pred), g.nodeFromId(l_k));
                            }
                            unsigned long l_k2 = l_k1 + 2 * sigma + 1;
                            edge = ae(g.nodeFromId(l_k1), g.nodeFromId(l_k2));
                            while (edge != INVALID && flow[edge] > 0.0) {
                                overflow = min(overflow, flow[edge]);
                                l_k1 = l_k2;
                                l_k2 += 2 * sigma + 1;
                                edge = ae(g.nodeFromId(l_k1), g.nodeFromId(l_k2));
                            }
                            flow[upper] -= overflow;
                            if (j1 < j2) {
                                unsigned long start = l_k;
                                unsigned long start_lower = l_k + 2 * (j2 - j1);
                                unsigned long u_start = l_k - 2 * j1 - 1;
                                SmartDigraph::Arc edge = ae(g.nodeFromId(u_start), g.nodeFromId(start));
                                flow[edge] -= overflow;
                                edge = ae(g.nodeFromId(u_start), g.nodeFromId(start_lower));
                                flow[edge] += overflow;
                                while (start != l_k) {
                                    unsigned long start_next = start + 2 * sigma + 1;
                                    unsigned long start_lower_next = start + 2 * sigma + 1;
                                    SmartDigraph::Arc edge = ae(g.nodeFromId(start), g.nodeFromId(start_next));
                                    flow[edge] -= overflow;
                                    edge = ae(g.nodeFromId(start_lower), g.nodeFromId(start_lower_next));
                                    flow[edge] += overflow;
                                    start = start_next;
                                    start_lower = start_lower_next;
                                }
                                unsigned long u_k_end = start - j1 * 2 - 1;
                                edge = ae(g.nodeFromId(start), g.nodeFromId(u_k_end));
                                flow[edge] -= overflow;
                                edge = ae(g.nodeFromId(start_lower), g.nodeFromId(start_lower + 2 * sigma + 1));
                                flow[edge] += overflow;
                                edge = ae(g.nodeFromId(u_k1), g.nodeFromId(start_lower + 2 * sigma + 1));
                                flow[edge] -= overflow;
                            } else {
                                unsigned long start = u_k + j2 * 2 + 1;
                                unsigned long start_lower = u_k + j1 * 2 + 1;
                                SmartDigraph::Arc edge = ae(g.nodeFromId(u_k), g.nodeFromId(start_lower));
                                flow[edge] -= overflow;
                                edge = ae(g.nodeFromId(start_lower), g.nodeFromId(start_lower + 2 * sigma + 1));
                                flow[edge] += overflow;
                                start += 2 * sigma + 1;
                                start_lower += 2 * sigma + 1;
                                while (start != l_k1) {
                                    unsigned long start_next = start + 2 * sigma + 1;
                                    unsigned long start_lower_next = start + 2 * sigma + 1;
                                    SmartDigraph::Arc edge = ae(g.nodeFromId(start), g.nodeFromId(start_next));
                                    flow[edge] -= overflow;
                                    edge = ae(g.nodeFromId(start_lower), g.nodeFromId(start_lower_next));
                                    flow[edge] += overflow;
                                    start = start_next;
                                    start_lower = start_lower_next;
                                }
                                unsigned long u_k_end = start - j2 * 2 - 1;
                                edge = ae(g.nodeFromId(start), g.nodeFromId(u_k_end));
                                flow[edge] -= overflow;
                                edge = ae(g.nodeFromId(start_lower), g.nodeFromId(u_k_end));
                                flow[edge] += overflow;
                            }
                        }
                    }
                }
                //round flow like increasing edges
                while (flow[upper] > (double) servers[server_counter].amount_servers) {
                    unsigned long l_k = u_k + 2 * sigma - 1;
                    unsigned long j = sigma;
                    for (int i = 0; i != sigma - 2; ++i) {
                        SmartDigraph::Arc edge = ae(g.nodeFromId(l_k), g.nodeFromId(u_k));
                        if (flow[edge] > 0.0) {
                            break;
                        } else {
                            j--;
                            l_k -= 2;
                        }
                    }
                    unsigned long l_ak = l_k + 1;
                    unsigned long l_k1 = l_k + 2 * sigma + 1;
                    SmartDigraph::Arc power_up = ae(g.nodeFromId(l_k), g.nodeFromId(u_k));
                    double overflow = min(flow[upper] - (double) servers[server_counter].amount_servers,
                                          flow[power_up]);
                    flow[upper] -= overflow;
                    //route flow to lower path
                    flow[power_up] -= overflow;
                    //route flow over lower paths
                    SmartDigraph::Arc lower = ae(g.nodeFromId(l_k), g.nodeFromId(l_ak));
                    SmartDigraph::Arc lower1 = ae(g.nodeFromId(l_ak), g.nodeFromId(l_k1));
                    SmartDigraph::Arc power_up1 = ae(g.nodeFromId(l_k1), g.nodeFromId(u_k1));
                    flow[lower] += overflow;
                    flow[lower1] += overflow;
                    flow[power_up1] += overflow;
                }
            }
            outer:;
            vertices_counter += 1 + 2 * sigma;
        }
    }
}

void pack_flow(ArcLookUp<SmartDigraph> &ae, const vector<Server> &servers, const unsigned long amnt_timesteps,
               const SmartDigraph &g,
               SmartDigraph::ArcMap<unsigned long> &cost, SmartDigraph::ArcMap<double> &flow) {
    //FIXME Make sure that no edges on the upper path are overwritten
    SmartDigraph::ArcMap<double> aux_flow(g);
    unsigned long server_counter = 0;
    unsigned long vertices_counter = 2 * amnt_timesteps;
    bool is_decreasing = false;
    double prior_flow = 0;
    unsigned long valley_start = 0;
    unsigned long valley_end = 0;
    unsigned long first_server = vertices_counter;
    while (server_counter < servers.size()) {
        unsigned long sigma = servers[server_counter].transition_costs.size();
        unsigned long u_k = vertices_counter;
        vertices_counter += 1 + 2 * sigma;
        unsigned long u_k1 = vertices_counter;
        SmartDigraph::Arc upper = ae(g.nodeFromId(u_k), g.nodeFromId((u_k1)));
        if (upper == INVALID) {
            vertices_counter += 1 + sigma;
            server_counter++;
            first_server = vertices_counter;
            continue;
        }
        if (flow[upper] < prior_flow) {
            is_decreasing = true;
            valley_start = u_k;
        } else if (flow[upper] > prior_flow && is_decreasing) {
            //matching pair found
            valley_end = u_k;
            double delta_k = 0;
            unsigned long l_k = valley_start + 1;
            vector<unsigned long> lower_flow;
            for (int i = 0; i != sigma; ++i) {
                SmartDigraph::Arc power_down = ae(g.nodeFromId(valley_start), g.nodeFromId((l_k)));
                if (flow[power_down] > 0) {
                    lower_flow.emplace_back(l_k);
                    delta_k += flow[power_down];
                }
                l_k += 2;
            }
            double delta_k_end = 0;
            unsigned long l_k_end = valley_end + 1;
            vector<unsigned long> lower_flow_end;
            for (int i = 0; i != sigma; ++i) {
                SmartDigraph::Arc power_up = ae(g.nodeFromId(valley_end), g.nodeFromId((l_k_end)));
                if (flow[power_up] > 0) {
                    delta_k_end += flow[power_up];
                    lower_flow_end.emplace_back((l_k_end));
                }
                l_k_end += 2;
            }
            double delta;
            vector<unsigned long> lower_paths;
            if (delta_k <= delta_k_end) {
                delta = delta_k;
                lower_paths = lower_flow;
            } else {
                delta = delta_k_end;
                lower_paths = lower_flow_end;
            }
            unsigned long cheapest_path = l_k;
            double min_cost = numeric_limits<double>::max();
            for (const unsigned long &j : lower_paths) {
                l_k = j;
                SmartDigraph::Arc power_down = ae(g.nodeFromId(valley_start), g.nodeFromId(l_k));
                double delta_j = flow[power_down];
                flow[power_down] = 0.0;

                unsigned long sum = 0;
                while (l_k < valley_end + 1 + 2 * sigma) {
                    unsigned long l_k1 = l_k + 1 + 2 * sigma;
                    SmartDigraph::Arc lower = ae(g.nodeFromId(l_k), g.nodeFromId(l_k1));
                    sum += cost[lower];
                    flow[lower] -= delta_j;
                    l_k = l_k1;
                }
                if (sum < min_cost) {
                    min_cost = sum;
                    cheapest_path = j;
                }
                SmartDigraph::Arc power_up = ae(g.nodeFromId(l_k), g.nodeFromId(valley_end));
                flow[power_up] -= delta_j;
            }
            l_k = cheapest_path;

            SmartDigraph::Arc power_down = ae(g.nodeFromId(valley_start), g.nodeFromId(l_k));
            aux_flow[power_down] += delta;
            while (l_k < valley_end + 1 + 2 * sigma) {
                unsigned long l_k1 = l_k + 1 + 2 * sigma;
                SmartDigraph::Arc lower = ae(g.nodeFromId(l_k), g.nodeFromId(l_k1));
                aux_flow[lower] += delta;
                l_k = l_k1;
            }
            SmartDigraph::Arc power_up = ae(g.nodeFromId(l_k), g.nodeFromId(valley_end));
            aux_flow[power_up] += delta;
            vertices_counter = first_server;
        }
    }
    for(SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        flow[a] = max(flow[a], aux_flow[a]);
    }
}

double mcmcf(bool is_debug, const SmartDigraph &g, const SmartDigraph::ArcMap<unsigned long> &capacity,
             const SmartDigraph::ArcMap<unsigned long> &capacity1,
             const SmartDigraph::ArcMap<unsigned long> &capacity2, SmartDigraph::ArcMap<unsigned long> &cost,
             SmartDigraph::NodeMap<long> &imbalances1,
             SmartDigraph::NodeMap<long> &imbalances2, Lp &lp, SmartDigraph::ArcMap<Lp::Col> &f1) {
    // Create an instance of the default LP solver
    // Add a column to the problem for each arc
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
        for (SmartDigraph::OutArcIt a(g, n); a != INVALID; ++a) {
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
    switch (lp.primalType()) {
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

    if (is_debug) {
        for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
            printf("%d: %.17g %.17g\n", g.id(a), (double) lp.primal(f1[a]), (double) lp.primal(f2[a]));
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
            for (auto it = buffer.begin(); it != buffer.end(); ++it) {
                if (!(*it).empty()) {
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



void benchmark(bool is_debug, int filenumber, uint64_t &generate, uint64_t &time_flow, string output) {
    srand(time(NULL));
    vector<Server> servers;
    vector<unsigned long> demands;
    read_test_file("tests/test_" + to_string(filenumber), servers, demands);
    ofstream file;
    ofstream file_rounding;
    ofstream data_file_rounding;
    ofstream data_file;
    file.open(output, std::ios_base::app);
    data_file.open(output + "_data", std::ios_base::app);
    file_rounding.open(output + "_rounding", std::ios_base::app);
    data_file_rounding.open(output + "_data_rounding", std::ios_base::app);
    long amount_nodes = 2 * (demands.size() + 1);
    long amount_edges = 0;

    for (auto it_servers = servers.begin(); it_servers != servers.end(); ++it_servers) {
        amount_nodes += (1 + 2 * it_servers->transition_costs.size()) * demands.size() + 1 +
                        it_servers->transition_costs.size();
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
    SmartDigraph::ArcMap<double> flow(g);

    Lp lp;
    SmartDigraph::ArcMap<Lp::Col> f1(g);

    auto start = chrono::high_resolution_clock::now();
    generate_graph(g, imbalances1, imbalances2, cost, capacity, capacity1, capacity2, servers, demands);
    auto start_flow = chrono::high_resolution_clock::now();

    double min_cost = mcmcf(is_debug, g, capacity, capacity1, capacity2, cost, imbalances1, imbalances2, lp, f1);
    auto end_min_flow = chrono::high_resolution_clock::now();

    ArcLookUp<SmartDigraph> ae(g);
    //ae.refresh();
    cout << "Starting Scaling" << std::endl;
    scale_flow(servers.size(), g, capacity1, flow, imbalances1, f1, lp);
    auto end_scale = chrono::high_resolution_clock::now();
    cout << "Starting Valleys" << std::endl;
    round_flow_valley(ae, servers, demands.size(), g, flow);
    auto end_valley = chrono::high_resolution_clock::now();
    cout << "Starting Increasing" << std::endl;
    round_flow_inc(ae, servers, demands.size(), g, flow);
    auto end_inc = chrono::high_resolution_clock::now();
    cout << "Starting Decreasing" << std::endl;
    round_flow_dec(ae, servers, demands.size(), g, flow);
    auto end_dec = chrono::high_resolution_clock::now();
    cout << "Starting peaks" << std::endl;
    round_peaks(ae,servers, demands.size(), g, flow);
    auto end_peaks = chrono::high_resolution_clock::now();
    cout << "Starting Reducing" << std::endl;
    reduce_to_m(ae, servers, demands.size(), g, capacity1, flow);
    auto end_reduce = chrono::high_resolution_clock::now();
    cout << "Starting Packing" << std::endl;
    pack_flow(ae, servers, demands.size(), g, cost, flow);
    auto end_pack = chrono::high_resolution_clock::now();

    cout << "Minimal Cost: " << min_cost << std::endl;

    auto scale_bench = chrono::duration_cast<std::chrono::nanoseconds>(end_scale-end_min_flow).count();
    auto valley_bench = chrono::duration_cast<std::chrono::nanoseconds>(end_valley-end_scale).count();
    auto inc_bench = chrono::duration_cast<std::chrono::nanoseconds>(end_inc-end_valley).count();
    auto dec_bench = chrono::duration_cast<std::chrono::nanoseconds>(end_dec-end_inc).count();
    auto peak_bench = chrono::duration_cast<std::chrono::nanoseconds>(end_peaks-end_dec).count();
    auto reduce_bench = chrono::duration_cast<std::chrono::nanoseconds>(end_reduce-end_peaks).count();
    auto pack_bench = chrono::duration_cast<std::chrono::nanoseconds>(end_pack-end_reduce).count();
    auto generate_bench = chrono::duration_cast<std::chrono::nanoseconds>(start_flow-start).count();
    auto flow_bench = chrono::duration_cast<std::chrono::nanoseconds>(end_min_flow-start_flow).count();
    uint64_t bench = generate_bench + flow_bench;
    generate += generate_bench;
    time_flow += flow_bench;
    file << (filenumber + 1) << " & " << amount_nodes << " & " << amount_edges << " & ";
    file << "\\SI{" << generate_bench << "}{\\nano\\second} & " << "\\SI{" << flow_bench << "}{\\nano\\second} & "
         << "\\SI{" << bench << "}{\\nano\\second}\\\\" << std::endl;
    file << "\\hline" << std::endl;
    file_rounding << (filenumber + 1) << " & " << amount_nodes << " & " << amount_edges << " & ";
    file_rounding << "\\SI{" << scale_bench << "}{\\nano\\second} & " << "\\SI{" << valley_bench << "}{\\nano\\second} & "
         << "\\SI{" << inc_bench << "}{\\nano\\second} &" << "\\SI{" << dec_bench << "}{\\nano\\second} & " << "\\SI{" << peak_bench << "}{\\nano\\second} & ";
    file_rounding << "\\SI{" << reduce_bench << "}{\\nano\\second} & " << "\\SI{" << pack_bench << "}{\\nano\\second}\\\\" << std::endl;
    file_rounding << "\\hline" << std::endl;
    data_file << generate_bench << " & " << flow_bench << " & " << bench << "\\\\" << std::endl;
    data_file_rounding << scale_bench << " & " << valley_bench << " & " << inc_bench << " & " << dec_bench << " & " << peak_bench << " & " << reduce_bench << " & " << pack_bench << "\\\\" << std::endl;

    file.close();
    file_rounding.close();
    data_file.close();
    data_file_rounding.close();
}

int main(int argc, char *argv[]) {
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
    for (int i = 0; i != amount_tests; ++i) {
        printf("running test %d\n", i);
        benchmark(false, i, generate, flow, "results/result");
    }

    printf("Ran %d tests.\n", amount_tests);
    printf("Overall time required for generating the graph: %" PRIu64 " nanoseconds\n", generate);
    printf("Overall time required for executing the mcmcf algorithm: %" PRIu64 " nanoseconds\n", flow);
    printf("Overall time required: %" PRIu64 " nanoseconds\n", generate + flow);
    return 0;
}