/// Problem Generator
/// 
/// 

use super::{grid, coords};

/// max(cos(x),0)
fn pcos(x:f64)->f64{
    let cos = x.cos();
    if cos > 0.0 { cos }
    else {0.0}
}

pub fn incident_flux(p: &grid::GridCell) -> f64 {
    let area = p.polygon.area();
    let centroid = p.polygon.center();
    centroid.theta.sin() * area * pcos(centroid.phi)
}

pub fn thermal_flux(p: &grid::GridCell) -> f64 {
    let area = p.polygon.area();
    let temperature = p.value;
    area * temperature.powi(4)
}

fn adv_flux_across_edge(edge: &coords::Edge,a: &grid::GridCell,b: &grid::GridCell) -> Result<f64,&'static str> {
    if !a.polygon.to_edges().contains(&edge) {
        Err("Edge not in polygon a")
    }
    else if !b.polygon.to_edges().contains(&edge) {
        Err("Edge not in polygon b")
    }
    else {
    Ok((b.value - a.value) * edge.phihat_dot_nhat())
    }
}

pub fn advective_flux(p: &grid::GridCell, network: &grid::GridNetwork) -> Result<f64,&'static str> {
    let mut flux = 0.0;
    let sides = p.polygon.to_edges();
    for side in sides.iter() {
        let neighbor_candidates = network.query_edge(side);
        if neighbor_candidates.len() > 2 { return Err("This edge separates more than two cells!");}
        let neighbor = {
            if neighbor_candidates[0] == p { neighbor_candidates[1] }
            else { neighbor_candidates[0] }
        };
        flux += adv_flux_across_edge(side,p,neighbor).unwrap();
    }
    Ok(flux)
}

fn diff_flux_across_edge(edge: &coords::Edge,a: &grid::GridCell,b: &grid::GridCell) -> Result<f64,&'static str> {
    let dist = a.polygon.center().angle_between(&b.polygon.center());
    let len_boundary = edge.len();
    Ok((b.value - a.value) / dist * len_boundary)
}

pub fn diffusive_flux(p: &grid::GridCell, network: &grid::GridNetwork) -> Result<f64,&'static str> {
    let mut flux = 0.0;
    let sides = p.polygon.to_edges();
    for side in sides.iter() {
        let neighbor_candidates = network.query_edge(side);
        if neighbor_candidates.len() > 2 { return Err("This edge separates more than two cells!");}
        let neighbor = {
            if neighbor_candidates[0] == p { neighbor_candidates[1] }
            else { neighbor_candidates[0] }
        };
        flux += diff_flux_across_edge(side,p,neighbor).unwrap();
    }
    Ok(flux)
}