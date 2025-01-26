use std::f64::consts::PI;


pub mod coords;
pub mod grid;
pub mod geometry;
pub mod pgen;
pub mod meshgen;

fn main() {
    let c1 = coords::Coordinate::new(0.0 ,PI/2.).unwrap();
    let c2 = coords::Coordinate::new(PI/16. ,PI/2.).unwrap();
    let c3 = coords::Coordinate::new(0.0 ,0.0).unwrap();
    let c4 = coords::Coordinate::new(PI/8.0 ,PI/2.).unwrap();
    let p1 = coords::Polygon::new(vec![c3.clone(),c1.clone(),c2.clone()]);
    let p2 = coords::Polygon::new(vec![c3.clone(),c2.clone(),c4.clone()]);
    // let g1 = grid::GridCell::new(p1.clone(),0.);
    // let g2 = grid::GridCell::new(p2.clone(),0.);
    // let e1 = coords::Edge::new(c2.clone(),c4.clone());
    // let net = grid::GridNetwork::new(vec![g1.clone(),g2.clone()]);
    // print!("{}",net.query_neighbors(g1)[0]==g2);
    print!("{}",p1.center());
}
