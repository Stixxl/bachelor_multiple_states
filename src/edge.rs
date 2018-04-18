pub struct Edge {
    pub cost: usize,
    pub capacity: f64,
    pub capacity_1: f64,
    pub flow_1: f64,
    pub capacity_2: f64,
    pub flow_2: f64
}

impl Edge {
    pub fn new(cost: usize, capacity: usize, capacity_1: usize, capacity_2: usize) -> Edge {
        Edge {
            cost,
            capacity: capacity as f64,
            capacity_1: capacity_1 as f64,
            capacity_2: capacity_2 as f64,
            flow_1 : 0.0,
            flow_2 : 0.0
        }
    }
}
