
use std::f64::consts::PI;
use std::iter::zip;
use core::fmt::Display;

use super::geometry::GreatCircle;

#[derive(PartialEq,Clone,Copy)]
pub struct Coordinate{
    pub phi:f64,
    pub theta:f64
}

impl Coordinate{
    pub fn new(phi:f64,theta:f64)->Result<Coordinate,&'static str>{
        if phi.is_nan() { Err("phi is nan") }
        else if theta.is_nan() { Err("theta is nan") }
        else if theta > PI { Err("theta must be less than PI") }
        else if theta < 0.0 { Err("theta must be greater than 0.0") }
        else { Ok(Coordinate{phi:phi,theta:theta}) }
    }
    pub fn from_cart(x:f64,y:f64,z:f64)->Result<Coordinate,&'static str>{
        if x.is_nan() { Err("x is nan in from_cart") }
        else if y.is_nan() { Err("y is nan in from_cart") }
        else if z.is_nan() { Err("z is nan in from_cart") }

        else { Ok(Coordinate::new((x).atan2(y),(z).acos()).unwrap()) }
    }
    pub fn cart(self:&Coordinate)-> Result<(f64,f64,f64),String> {
        let x = self.phi.cos()*self.theta.sin();
        let y = self.phi.sin()*self.theta.sin();
        let z = self.theta.cos();
        if x.is_nan() { Err(String::from("x is nan in cart ") + self.to_string().as_str()) }
        else if y.is_nan() { Err(String::from("x is nan in cart ") + self.to_string().as_str()) }
        else if z.is_nan() { Err(String::from("x is nan in cart ") + self.to_string().as_str()) }
        else { Ok((x,y,z)) }
    }

    pub fn dot(self:&Coordinate,other:&Coordinate)->Result<f64,&'static str> {
        let res = self.theta.sin() * other.theta.sin() * (self.phi - other.phi).cos() + self.theta.cos()*other.theta.cos();
        if res.is_nan() { Err("NaN in dot product")}
        else { Ok(res) }
    }
    pub fn cross(self:&Coordinate,other:&Coordinate)->Result<Coordinate,&'static str> {
        let (x1,y1,z1) = self.cart().unwrap();
        let (x2,y2,z2) = other.cart().unwrap();
        let x = y1*z2 - y2*z1;
        let y = x2*z1 - x1*z2;
        let z = x1*y2 - x2*y1;
        Coordinate::from_cart(x,y,z)
    }

    pub fn sin_dist(self:&Coordinate,other:&Coordinate)->f64 {
        let cross = self.cross(other).unwrap();
        let (x,y,z) = cross.cart().unwrap();
        (x*x + y*y + z*z).sqrt()
    }
    pub fn angle_between(self:&Coordinate,other:&Coordinate)->f64 {
        if self.dot(other).unwrap() < 0.0 {
            let (x1,y1,z1) = self.cart().unwrap();
            let (x2,y2,z2) = other.cart().unwrap();
            let dx = -x1 - x2;
            let dy = -y1 - y2;
            let dz = -z1 - z2;
            let half_mag = (dx*dx + dy*dy + dz*dz).sqrt()/2.0;
            
            return PI - 2.0 * half_mag.asin()
        }
        else {
            let (x1,y1,z1) = self.cart().unwrap();
            let (x2,y2,z2) = other.cart().unwrap();
            let dx = x1 - x2;
            let dy = y1 - y2;
            let dz = z1 - z2;
            let half_mag = (dx*dx + dy*dy + dz*dz).sqrt()/2.0;

            return 2.0 * half_mag.asin()
        }
    }
}

impl Display for Coordinate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.phi, self.theta)
    }
}

pub fn spherical_angle(
    A: &Coordinate,
    B: &Coordinate,
    C: &Coordinate
) -> f64 {
    // Measure the angle of AB and BC
    
    // b is the distance AC
    let cos_b = A.dot(C).unwrap();

    // c is the distance AB
    let cos_c = A.dot(B).unwrap();
    let sin_c = A.sin_dist(B);

    // a is the distance BC
    let cos_a = B.dot(C).unwrap();
    let sin_a = B.sin_dist(C);

    let cos_B = (cos_b - cos_c * cos_a)/(sin_c*sin_a);
    cos_B.acos()
}

pub fn midpoint(a: &Coordinate, b: &Coordinate) -> Result<Coordinate, &'static str> {
    let (x1,y1,z1) = a.cart().unwrap();
    let (x2,y2,z2) = b.cart().unwrap();
    let x = (x1 + x2)/2.0;
    let y = (y1 + y2)/2.0;
    let z = (z1 + z2)/2.0;
    Coordinate::from_cart(x,y,z)
}

pub fn length(a: &Coordinate, b: &Coordinate) -> f64 {
    a.angle_between(b)
}
pub fn phihat_dot_nhat(a: &Coordinate, nhat: &Coordinate) -> f64 {
    let (nx, ny, _) = nhat.cart().unwrap();
    let phi = a.phi;
    ny * phi.cos() - nx * phi.sin()
}

pub struct Edge{
    pub a:Coordinate,
    pub b:Coordinate
}

impl Edge{
    pub fn new(a:Coordinate,b:Coordinate)->Edge{
        Edge{a:a,b:b}
    }
    pub fn len(self:&Edge)->f64{
        length(&self.a,&self.b)
    }
    pub fn midpoint(self:&Edge)->Result<Coordinate,&'static str>{
        midpoint(&self.a,&self.b)
    }
    pub fn get_great_circle(self:&Edge)->GreatCircle{
        GreatCircle::from_coords(self.a.clone(),self.b.clone())
    }
    /// Compute $\hat{\phi} \cdot \hat{n}$ using Simpson's rule
    pub fn phihat_dot_nhat(self:&Edge)->f64{
        let nhat = self.get_great_circle().nhat();
        let f_a = phihat_dot_nhat(&self.a,&nhat);
        let f_b = phihat_dot_nhat(&self.b,&nhat);
        let f_mid = phihat_dot_nhat(&self.midpoint().unwrap(),&nhat);
        self.len()/6.0 * (f_a + f_b + 4.0 * f_mid)
    }
}
impl PartialEq for Edge{
    fn eq(self:&Edge,other:&Edge)->bool{
        (self.a == other.a && self.b == other.b) || (self.a == other.b && self.b == other.a)
    }
}

pub struct PolyLine{
    nodes:Vec<Coordinate>
}

impl PolyLine{
    pub fn new(nodes:Vec<Coordinate>)->PolyLine{
        PolyLine{nodes:nodes}
    }
    pub fn to_edges(self:&PolyLine)->Vec<Edge>{
        let mut edges:Vec<Edge> = Vec::new();
        for i in 0..self.nodes.len()-1{
            edges.push(Edge::new(self.nodes[i].clone(),self.nodes[i+1].clone()));
        }
        edges
    }
}
#[derive(Clone)]
pub struct Polygon{
    pub nodes:Vec<Coordinate>
}

impl Polygon{
    pub fn new(nodes:Vec<Coordinate>)->Polygon{
        Polygon{nodes:nodes}
    }
    pub fn to_edges(self:&Polygon)->Vec<Edge>{
        let mut edges:Vec<Edge> = Vec::new();
        for i in 0..self.nodes.len()-1{
            edges.push(Edge::new(self.nodes[i].clone(),self.nodes[i+1].clone()));
        }
        edges.push(Edge::new(self.nodes[self.nodes.len()-1].clone(),self.nodes[0].clone()));
        edges
    }
    pub fn interior_angles(self:&Polygon)->Vec<f64>{
        let mut angles:Vec<f64> = Vec::new();
        for i in 0..self.nodes.len()-2{
            angles.push(spherical_angle(&self.nodes[i],&self.nodes[i+1],&self.nodes[i+2]));
        }
        angles.push(spherical_angle(&self.nodes[self.nodes.len()-2],&self.nodes[self.nodes.len()-1],&self.nodes[0]));
        angles.push(spherical_angle(&self.nodes[self.nodes.len()-1],&self.nodes[0],&self.nodes[1]));
        angles
    }

    pub fn area(self:&Polygon)->f64{
        let mut area:f64 = 0.0;
        let angles = self.interior_angles();
        let n = angles.len();
        for i in 0..n{
            area += angles[i];
        }
        area - (n as f64 - 2.0) * PI 
    }
    pub fn center(self:&Polygon)->Coordinate{
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        for edge in self.to_edges(){
            let a = edge.a;
            let b = edge.b;
            let v = a.cross(&b).unwrap();
            let angle = a.angle_between(&b);
            let (_x, _y, _z) = v.cart().unwrap();
            x += _x * angle;
            y += _y * angle;
            z += _z * angle;
        }
        let mag = (x*x + y*y + z*z).sqrt();
        Coordinate::from_cart(x/mag,y/mag,z/mag).unwrap()
    }
}

impl Display for Polygon{
    fn fmt(self:&Polygon,f:&mut std::fmt::Formatter)->std::fmt::Result{
        let mut s = String::new();
        s += "Polygon: [";
        for node in self.nodes.iter(){
            s += &node.to_string();
            s += ", ";
        }
        s += "]";
        write!(f,"{}",s)
    }
}

impl PartialEq for Polygon{
    fn eq(self:&Polygon,other:&Polygon)->bool{
        if self.nodes.len() != other.nodes.len(){
            return false;
        }
        let node0 = self.nodes[0];
        let i0 :isize = {
            let mut found = false;
            let mut _i0 = 0;
            for i in 0..other.nodes.len(){
                if other.nodes[i] == node0{
                    found = true;
                    _i0 = i;
                    break;
                }
            }
            if found{
                _i0 as isize
            }
            else{
                -1
            }
        };
        if i0 == -1{
            return false;
        }
        let i0 = i0 as usize;
        for (i,j) in zip(0..self.nodes.len(), i0..other.nodes.len() +i0){
            let _j = j % other.nodes.len();
            if self.nodes[i] != other.nodes[_j]{
                return false;
            }
        }
        true
    }
}