use core::f64;
use std::{error, f64::consts::PI};

use log::{info,error};
use std::time::{Instant,Duration};

/// Problem Generator
/// 
/// 

use super::{grid, coords,meshgen::icoshedron};

static COURANT_NUMBER: f64 = 0.4;

pub enum CFL_Limiter {
    AdvectionLimited(f64),
    DiffusionLimited(f64),
    SourceLimited(f64),
    NoLimit(f64)
}

pub fn get_timestep(eps1:f64,eps2:f64,max_temp:f64,dx:f64)->CFL_Limiter {
    let adv = {
        if eps1 == 0.0 { f64::INFINITY }
        else {COURANT_NUMBER * dx / eps1}
    };
    let diff = {
        if eps2 == 0.0 { f64::INFINITY }
        else {COURANT_NUMBER * dx.powi(2) / eps2 / 2.0}
    };
    let source = {
        if max_temp == 0.0 { f64::INFINITY }
        else { COURANT_NUMBER / 4.0 / max_temp.powi(3)}
    };

    if adv < diff && adv < source { CFL_Limiter::AdvectionLimited(adv) }
    else if diff < adv && diff < source { CFL_Limiter::DiffusionLimited(diff) }
    else if source < adv && source < diff { CFL_Limiter::SourceLimited(source) }
    else { CFL_Limiter::NoLimit(1.0) }
}

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
        let nhat = edge.get_great_circle().nhat();
        let nhat = {
            let (ax,ay,az) = a.polygon.center().cart().unwrap();
            let (bx,by,bz) = b.polygon.center().cart().unwrap();
            let dx = bx - ax;
            let dy = by - ay;
            let dz = bz - az;
            let (nx,ny,nz) = nhat.cart().unwrap();
            let _dot = dx*nx + dy*ny + dz*nz;
            if _dot < 0.0 {
                coords::Coordinate::from_cart(-nx, -ny, -nz).unwrap()
            }
            else { nhat }
        };
        let upwind_i = {
            if coords::phihat_dot_nhat(&edge.a,&nhat) > 0.0 { 
                a
             } else { 
                b
              }
        };
        let upwind_m = {
            if coords::phihat_dot_nhat(&edge.midpoint().unwrap(),&nhat) > 0.0 { a } else { b }
        };
        let upwind_f = {
            if coords::phihat_dot_nhat(&edge.b,&nhat) > 0.0 { a } else { b }
        };

        let f_i = upwind_i.value * edge.a.theta.sin() * coords::phihat_dot_nhat(&edge.a, &nhat);
        let f_m = upwind_m.value * edge.midpoint().unwrap().theta.sin() * coords::phihat_dot_nhat(&edge.midpoint().unwrap(), &nhat);
        let f_f = upwind_f.value * edge.b.theta.sin() * coords::phihat_dot_nhat(&edge.b, &nhat);
        
        Ok(edge.len() / 6.0 * (f_i + 4.0*f_m + f_f))  
    }
}

pub fn advective_flux(p: &grid::GridCell, network: &grid::GridNetwork) -> Result<f64,&'static str> {
    let mut flux = 0.0;
    let sides = p.polygon.to_edges();
    for side in sides.iter() {
        let neighbor_candidates = network.query_edge(side);
        if neighbor_candidates.len() > 2 { 
            let mut s = String::from("This edge separates more than two cells!");
            s += &format!("\nEdge from {:}",side.a);
            s += &format!("\n to       {:}",side.b);
            s += &format!("\nNeighbors:");
            for can in neighbor_candidates.iter() {
                s += &format!("\n{:}",can.polygon);
            }
            error!("{}",s);
            return Err("This edge separates more than two cells!");
        }
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
        if neighbor_candidates.len() > 2 { 
            let mut s = String::from("This edge separates more than two cells!");
            s += &format!("\nEdge from {:}",side.a);
            s += &format!("\n to       {:}",side.b);
            s += &format!("\nNeighbors:");
            for can in neighbor_candidates.iter() {
                s += &format!("\n{:}",can.polygon);
            }
            error!("{}",s);
            return Err("This edge separates more than two cells!");
        }
        let neighbor = {
            if neighbor_candidates[0] == p { neighbor_candidates[1] }
            else { neighbor_candidates[0] }
        };
        flux += diff_flux_across_edge(side,p,neighbor).unwrap();
    }
    Ok(flux)
}

fn format_debug_output(p: &grid::GridCell, network: &grid::GridNetwork) -> String {
    let mut s = String::from("Cell\n");
    s += "====\n";
    s += &format!("Value: {}\n",p.value);
    s += &format!("Area:  {}\n",p.polygon.area());
    s += &format!("Vertices:\n");
    for v in p.polygon.nodes.iter() {
        s += &format!("\t{:2}\n",v);
    }
    let neighbors = network.query_neighbors(&p);
    s += &format!("Neighbors:\n");
    for (i, n) in neighbors.iter().enumerate() {
        s += &format!("{})\n",i);
        s += &format!("Value: {}\n",n.value);
        s += &format!("Area:  {}\n",n.polygon.area());
        s += &format!("Vertices:\n");
        for v in n.polygon.nodes.iter() {
            s += &format!("\t{:2}\n",v);
        }
    }
    s
}



pub fn get_next_value(p: &grid::GridCell, network: &grid::GridNetwork,eps1: f64, eps2: f64, dt: f64) -> Result<f64,&'static str> {
    let area = p.polygon.area();
    let _incident_flux = incident_flux(p) * dt / area;
    let _thermal_flux = -thermal_flux(p) * dt / area;
    let _advective_flux = -advective_flux(p,network)? * dt / area * eps1;
    let _diffusive_flux = diffusive_flux(p,network)? * dt / area * eps2;
    let next_value = p.value + _incident_flux + _thermal_flux + _advective_flux + _diffusive_flux;
    if next_value < 0.0 {
        let mut s = String::from("Negative temperature in cell update");
        s += &format!("\n{:}\n\n",format_debug_output(p,network));
        s += &format!("\nIncident flux: {}",_incident_flux);
        s += &format!("\nThermal flux: {}",_thermal_flux);
        s += &format!("\nAdvective flux: {}",_advective_flux);
        s += &format!("\nDiffusive flux: {}",_diffusive_flux);
        s += &format!("\nNext value: {}",next_value);
        error!("{}",s);

        return Err("Negative temperature");
    }
    Ok(next_value)
}

pub enum InitialCondition {
    Constant(f64),
    Radiative
}

pub fn init_mesh(subdivisions: u32, initial_condition: InitialCondition) -> grid::GridNetwork {
    let polygons = icoshedron(subdivisions);
    let mut cells: Vec<grid::GridCell> = Vec::new();
    for p in polygons.iter() {
        let value = match initial_condition {
            InitialCondition::Constant(c) => c,
            InitialCondition::Radiative => pcos(p.center().phi) * p.center().theta.sin()
        };
        cells.push(grid::GridCell::new(p.clone(),value));
    }
    grid::GridNetwork::new(cells)
}

pub fn get_next_mesh(network: grid::GridNetwork, eps1: f64, eps2: f64) -> grid::GridNetwork {
    let max_temp = network.max_value();
    let dx = (4.0 * PI / network.cells.len() as f64).sqrt();
    let dt_result = get_timestep(eps1, eps2, max_temp, dx);
    let dt =match dt_result {
        CFL_Limiter::AdvectionLimited(adv) => {info!("Advection limited timestep: {}",adv); adv},
        CFL_Limiter::DiffusionLimited(diff) => {info!("Diffusion limited timestep: {}",diff); diff},
        CFL_Limiter::SourceLimited(source) => {info!("Source limited timestep: {}",source);source},
        CFL_Limiter::NoLimit(no_limit) => {info!("No limit timestep: {}",no_limit);no_limit}
    };
    let mut new_cells: Vec<grid::GridCell> = Vec::new();
    let start_time = Instant::now();
    for cell in network.cells.iter() {
        let value = get_next_value(cell,&network,eps1,eps2,dt).unwrap();
        new_cells.push(grid::GridCell::new(cell.polygon.clone(),value));
    }
    let end_time = Instant::now();
    info!("Mesh update time: {:?}",end_time - start_time);
    grid::GridNetwork::new(new_cells)

}

pub enum EnergyBalance {
    Balanced(f64),
    ExcessIncident(f64),
    ExcessThermal(f64),
}

pub fn check_energy_balance(network: &grid::GridNetwork) -> EnergyBalance {
    let mut energy_in = 0.0;
    let mut energy_out = 0.0;
    for cell in network.cells.iter() {
        energy_in += incident_flux(cell);
        energy_out += thermal_flux(cell);
    }
    let excess = energy_in - energy_out;
    let excess_as_pct = excess / energy_in * 100.0;
    if excess_as_pct.abs() < 1.0 {
        log::info!("Energy balance: {}%",excess_as_pct);
        EnergyBalance::Balanced(excess_as_pct)
    }
    else if excess_as_pct > 0.0 {
        log::warn!("Excess incident energy: {}%",excess_as_pct);
        EnergyBalance::ExcessIncident(excess_as_pct)
    }
    else {
        log::warn!("Excess thermal energy: {}%",excess_as_pct);
        EnergyBalance::ExcessThermal(excess_as_pct)
    }
}