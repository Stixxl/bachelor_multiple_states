mod server;
mod edge;

use server::Server;
use edge::Edge;

use std::collections::HashMap;

use std::fs::File;
use std::io::Write;

fn main() {
    let mut servers: Vec<Server> = Vec::new();
    let mut demands: Vec<usize> = Vec::new();

    let consumption_rate = vec![vec![3,3,3], vec![2,2,2], vec![1,1,1]];
    let transition_costs = vec![1,2];
    let server = Server {
        consumption_rate: consumption_rate,
        transition_costs: transition_costs,
        amount_servers: 3
    };
    servers.push(server);

    let consumption_rate = vec![vec![5,5,5], vec![3,3,3], vec![1,1,1]];
    let transition_costs = vec![2,4];
    let server = Server {
        consumption_rate: consumption_rate,
        transition_costs: transition_costs,
        amount_servers: 3
    };
    servers.push(server);

    demands.push(1);
    demands.push(2);
    let (graph, imbalances_1, imbalances_2) = construct_graph(servers, demands);
    println!("Imbalances for commodity 1: {:?}", imbalances_1);
    println!("Imbalances for commodity 2: {:?}", imbalances_2);
}

fn construct_graph(servers: Vec<Server>, demands: Vec<usize>) -> (HashMap<usize, HashMap<usize, Edge>>, Vec<isize>, Vec<isize>) {
    let mut graph: HashMap<usize, HashMap<usize, Edge>> = HashMap::new();
    let mut size = 2 * demands.len() + 1; //FIXME set to actual amount of vertices
    let mut vertices_counter: usize = 0;

    let mut amount_servers = 0;
    let mut d_k = 0;

    for server in servers.iter() {
        d_k += server.amount_servers * (server.transition_costs.len() - 1);
        amount_servers += server.amount_servers;
        size += (1 + 2 * server.transition_costs.len()) * demands.len() + 1 + server.transition_costs.len();
    }
    let mut imbalances_1: Vec<isize> = vec![0; size];
    let mut imbalances_2: Vec<isize> = vec![0; size];

    let a_0 = vertices_counter;
    imbalances_1[a_0] = amount_servers as isize;
    vertices_counter += 1;

    let b_0 = vertices_counter;
    imbalances_1[b_0] = -1 * amount_servers as isize;
    vertices_counter += 1;

    for demand in demands.iter() {
        let a_k = vertices_counter;
        imbalances_2[a_k] = (d_k + demand) as isize;
        vertices_counter += 1;

        let b_k = vertices_counter;
        imbalances_2[b_k] = -1 * (d_k + demand) as isize;
        vertices_counter += 1;
    }

    for server in servers.iter() {
        //TODO _1handle first timestep differently no edges from lower paths to upper paths except
        //for lowest; same for last
        for j in 0..demands.len() {
            let a_k = 2 * (j+1);
            let b_k = a_k + 1;

            let u_k = vertices_counter;
            vertices_counter += 1;
            for (k, transition_cost) in server.transition_costs.iter().enumerate() {

                let l_km = vertices_counter;
                vertices_counter += 1;
                let l_kma = vertices_counter;
                vertices_counter += 1;
                let mut l_km1 = vertices_counter + 2 * server.transition_costs.len() - 1;
                if j == demands.len() - 1 {
                    l_km1 -= k;
                }
                if j == 0 && k == server.transition_costs.len() - 1 {
                    graph.entry(a_0).or_insert(HashMap::new()).insert(l_km, Edge::new(0, server.amount_servers, server.amount_servers, 0));
                    graph.entry(l_km).or_insert(HashMap::new()).insert(u_k, Edge::new(*transition_cost, server.amount_servers, server.amount_servers, 0));
                }
                if j != 0 {
                    graph.entry(l_km).or_insert(HashMap::new()).insert(u_k, Edge::new(*transition_cost, server.amount_servers, server.amount_servers, 0));
                    graph.entry(u_k).or_insert(HashMap::new()).insert(l_km, Edge::new(0, server.amount_servers, server.amount_servers, 0));
                }
                graph.entry(l_km).or_insert(HashMap::new()).insert(l_kma, Edge::new(server.consumption_rate[k+1][j], server.amount_servers, server.amount_servers, 0));
                graph.entry(l_kma).or_insert(HashMap::new()).insert(l_km1, Edge::new(0, server.amount_servers, server.amount_servers, server.amount_servers));
                graph.entry(a_k).or_insert(HashMap::new()).insert(l_kma, Edge::new(0, server.amount_servers, 0, server.amount_servers));
                graph.entry(l_km1).or_insert(HashMap::new()).insert(b_k, Edge::new(0, server.amount_servers, 0, server.amount_servers));
                }
            let u_k1 = vertices_counter;
            graph.entry(u_k).or_insert(HashMap::new()).insert(u_k1, Edge::new(server.consumption_rate[0][j], server.amount_servers, server.amount_servers, 0));
        }
        //handle last timestep
        let u_kn = vertices_counter;
        //increase counter by all intermediate lower paths vertices l_{i,j,k}
        vertices_counter += server.transition_costs.len();
        let l_kn = vertices_counter;
        vertices_counter += 1;

        graph.entry(u_kn).or_insert(HashMap::new()).insert(l_kn, Edge::new(0, server.amount_servers, server.amount_servers, 0));
        graph.entry(l_kn).or_insert(HashMap::new()).insert(b_0, Edge::new(0, server.amount_servers, server.amount_servers, 0));

    }
    write_graph_in_DOT(&graph, &imbalances_1, &imbalances_2, demands.len());
    (graph, imbalances_1, imbalances_2)
}
fn write_graph_in_DOT(graph: &HashMap<usize, HashMap<usize, Edge>>, imbalances_1: &Vec<isize>, imbalances_2: &Vec<isize>, amount_timesteps: usize) {
    let mut dot_text = String::new();
    dot_text.push_str("digraph {\n");
    dot_text.push_str("rankdir=LR;\n");
    dot_text.push_str("node [shape=circle, fixedsize=true, width=0.25, color=\"black\",fillcolor=white, style=\"filled,solid\", fontcolor=red, fontsize=12]\nedge [ penwidth=0.75, color=black ]\n");

    for (node, edges) in graph.iter() {

        for (end_node, edge) in edges.iter() {
            dot_text.push_str(&node.to_string());
            dot_text.push_str(" -> ");
            dot_text.push_str(&end_node.to_string());
            if edge.cost > 0 {
                dot_text.push_str(" [ label=\"");
                dot_text.push_str(&edge.cost.to_string());
                dot_text.push_str("\" ]");
            }
            dot_text.push_str(";");
            dot_text.push_str("\n");
        }
    }

    dot_text.push_str("}");
    let mut file = File::create("/home/corbinian/Documents/University/WS17_18/Bachelorthesis/multiple_states/tmp/graph.gv").expect("Could not create file!");
    file.write_all(dot_text.as_bytes()).expect("Could not write data!");
}

