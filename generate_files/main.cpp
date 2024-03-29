#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <string>

using namespace std;
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

void generate_test(int testnumber, int amount_servers, int amount_demands, int max_transition, int max_consumption, int max_per_server, int max_lower) {
    ofstream file;
    file.open("tests/test_" + to_string(testnumber));
    srand(time(nullptr));
    int total_amount_servers = 0;

    vector<unsigned int> demands;
    vector<Server> servers;

    file << amount_servers << " " << amount_demands << std::endl;
    for(int i = 0; i != amount_servers; ++i) {
        int amount_per_server = (rand() % max_per_server) + 1;
        vector<unsigned long> tc;
        vector<vector<unsigned long>> cr;
        total_amount_servers += amount_per_server;
        int amount_lower = (rand() % max_lower) + 1;
        for (int k = 0; k != amount_lower + 1; ++k) {
            vector<unsigned long> crk;
            if(k != 0) {
                tc.push_back((rand() % max_transition) + 1);
            }
        for (int j = 0; j != amount_demands; ++j) {
                unsigned int consumption_rate;
            if (k == 0) {
                consumption_rate = (rand() % (max_consumption - amount_lower)) + amount_lower;
            } else {
                consumption_rate = (rand() % (cr[k-1][j] - amount_lower + k)) + amount_lower - k;
        }
        crk.push_back(consumption_rate);
        }
        cr.push_back(crk);
    }
    servers.push_back(Server(cr, tc, amount_per_server));
    }

    for(int i = 0; i != amount_demands; ++i) {
        int demand = rand() % (total_amount_servers + 1);
        demands.push_back(demand);
    }
    for (unsigned int &demand : demands) {
        file << demand << " ";
    }
    for (auto &server : servers) {
        file << std::endl;
        file << server.amount_servers << " " << server.transition_costs.size() << std::endl;
        for (unsigned long &transition_cost : server.transition_costs) {
         file << transition_cost << " ";
        }
        for (auto &it_consumption_rates : server.consumption_rate) {
            file << std::endl;
            for (unsigned long &it_consumption_rate : it_consumption_rates) {
                file << it_consumption_rate << " ";
            }
        }
    }

    printf("Creating an test instance with %d servers and %d timesteps; max transition cost %d; max consumption rate: %d; max servers per class: %d; maxi amount lower paths %d",
           amount_servers, amount_demands, max_transition, max_consumption, max_per_server, max_lower);
    file.close();

}

int main(int argc, char * argv[]) {
    string testnumber = argv[1];
    string amount_servers = argv[2];
    string amount_demands = argv[3];
    string max_transition = argv[4];
    string max_consumption = argv[5];
    string max_per_server = argv[6];
    string max_lower = argv[7];
    generate_test(stoi(testnumber), stoi(amount_servers), stoi(amount_demands), stoi(max_transition), stoi(max_consumption), stoi(max_per_server), stoi(max_lower));
    return 0;
}