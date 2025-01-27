use std::f64::consts::PI;

use log::info;
use simple_logger::{SimpleLogger};


pub mod coords;
pub mod grid;
pub mod geometry;
pub mod pgen;
pub mod meshgen;

fn main() {
    SimpleLogger::new().with_level(log::LevelFilter::Debug).init().unwrap();
    let n = 3;
    let e1 = 1.0;
    let e2 = 1.0;
    let mut net = pgen::init_mesh(n, pgen::InitialCondition::Constant(0.0));
    print!("Maximum value of the mesh is: {}", net.max_value());
    for i in 0..100 {
        net = pgen::get_next_mesh(net, e1, e2);
        _ = pgen::check_energy_balance(&net);
        info!("Maximum value of the mesh is: {}", net.max_value());
        info!("Average value of the mesh is: {}", net.average_value());
        info!("Minimum value of the mesh is: {}", net.min_value());
        
    }
}
