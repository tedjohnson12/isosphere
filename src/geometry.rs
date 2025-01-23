
use std::f32::consts::PI;

use super::coords;

pub struct GreatCircle {
    x: f64,
    y: f64,
    z: f64,
}

impl GreatCircle {
    pub fn new(x: f64, y: f64, z: f64) -> GreatCircle {
        let mag = (x*x + y*y + z*z).sqrt();
        
        GreatCircle{ x:x/mag, y:y/mag, z:z/mag }
    }
    pub fn from_coords(a: coords::Coordinate, b: coords::Coordinate) -> GreatCircle {
        let normal = a.cross(&b).unwrap();
        let mag = normal.dot(&normal).unwrap().sqrt();    
        let (x,y,z) = normal.cart().unwrap();
        GreatCircle{ x:x/mag, y:y/mag, z:z/mag }
    }
    pub fn nhat(&self) -> coords::Coordinate {
        coords::Coordinate::from_cart(self.x, self.y, self.z).unwrap()
    }
    pub fn inclination(&self) -> f64 {
        let zhat = coords::Coordinate::from_cart(0.0,0.0,1.0).unwrap();
        self.nhat().dot(&zhat).unwrap().acos()
    }

}
