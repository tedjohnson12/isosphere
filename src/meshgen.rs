//! Generate a triangular mesh
//! 
//! 

use std::{f64::consts::PI};

use crate::coords::subdivide_polygon;
use log::{info,warn};

use super::coords;


/// From https://danielsieger.com/blog/2021/01/03/generating-platonic-solids.html
pub fn icoshedron(n_subdivisions: u32) -> Vec<coords::Polygon> {
    info!("Generating icoshedron with {} subdivisions",n_subdivisions);
    let golden_ratio = (1. + (5. as f64).sqrt()) / 2.0;
    let a = 1.0;
    let b = 1.0/golden_ratio;
    let mag = (a*a + b*b).sqrt();
    let a = a/mag;
    let b = b/mag;

    let v1 = coords::Coordinate::from_cart(0.0,b,-a).unwrap();
    let v2 = coords::Coordinate::from_cart(b,a,0.0).unwrap();
    let v3 = coords::Coordinate::from_cart(-b,a,0.0).unwrap();
    let v4 = coords::Coordinate::from_cart(0.0,b,a).unwrap();
    let v5 = coords::Coordinate::from_cart(0.0,-b,a).unwrap();
    let v6 = coords::Coordinate::from_cart(-a,0.0,b).unwrap();
    let v7 = coords::Coordinate::from_cart(0.0,-b,-a).unwrap();
    let v8 = coords::Coordinate::from_cart(a,0.0,-b).unwrap();
    let v9 = coords::Coordinate::from_cart(a,0.0,b).unwrap();
    let v10 = coords::Coordinate::from_cart(-a,0.0,-b).unwrap();
    let v11 = coords::Coordinate::from_cart(b,-a,0.0).unwrap();
    let v12 = coords::Coordinate::from_cart(-b,-a,0.0).unwrap();

    let cells = vec![
        coords::Polygon::new(vec![v3.clone(),v2.clone(),v1.clone()]),
        coords::Polygon::new(vec![v2.clone(),v3.clone(),v4.clone()]),
        coords::Polygon::new(vec![v6.clone(),v5.clone(),v4.clone()]),
        coords::Polygon::new(vec![v5.clone(),v9.clone(),v4.clone()]),
        coords::Polygon::new(vec![v8.clone(),v7.clone(),v1.clone()]),
        coords::Polygon::new(vec![v7.clone(),v10.clone(),v1.clone()]),
        coords::Polygon::new(vec![v12.clone(),v11.clone(),v5.clone()]),
        coords::Polygon::new(vec![v11.clone(),v12.clone(),v7.clone()]),
        coords::Polygon::new(vec![v10.clone(),v6.clone(),v3.clone()]),
        coords::Polygon::new(vec![v6.clone(),v10.clone(),v12.clone()]),
        coords::Polygon::new(vec![v9.clone(),v8.clone(),v2.clone()]),
        coords::Polygon::new(vec![v8.clone(),v9.clone(),v11.clone()]),
        coords::Polygon::new(vec![v3.clone(),v6.clone(),v4.clone()]),
        coords::Polygon::new(vec![v9.clone(),v2.clone(),v4.clone()]),
        coords::Polygon::new(vec![v10.clone(),v3.clone(),v1.clone()]),
        coords::Polygon::new(vec![v2.clone(),v8.clone(),v1.clone()]),
        coords::Polygon::new(vec![v12.clone(),v10.clone(),v7.clone()]),
        coords::Polygon::new(vec![v8.clone(),v11.clone(),v7.clone()]),
        coords::Polygon::new(vec![v6.clone(),v12.clone(),v5.clone()]),
        coords::Polygon::new(vec![v11.clone(),v9.clone(),v5.clone()]),
    ];
    if n_subdivisions == 0 {
        cells
    }
    else {
        let mut cells = cells;
        let n_cells = cells.len();
        let expected_area = 4.0 * PI / (n_cells as f64)/3.0;
        for i in 0..n_subdivisions {
            info!("Starting subdivision {}",i);
            let mut new_cells:Vec<coords::Polygon> = Vec::new();
            for cell in cells.iter() {
                let subdivisions = subdivide_polygon(cell.clone());
                for s in subdivisions.iter() {
                    let area = s.area();
                    if 1.0 - area / expected_area > 0.01 {
                        warn!("Subdivided cell area is off by {}%",(1.0 - area / expected_area) * 100.0);
                    }



                    new_cells.push(s.clone());
                }
            }
            cells = new_cells;
            info!("There are now {} cells",cells.len());
        }
        cells
    }

}