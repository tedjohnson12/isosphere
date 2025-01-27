
use std::fmt::Display;
use super::coords;

#[derive(Clone)]
pub struct GridCell{
    pub polygon:coords::Polygon,
    pub value:f64
}

impl GridCell{
    pub fn new(polygon:coords::Polygon,value:f64)->GridCell{
        if value < 0.0 {panic!{"Cell value cannot be negative"}}
        GridCell{polygon:polygon,value:value}
    }
}

impl PartialEq for GridCell{
    fn eq(self:&GridCell,other:&GridCell)->bool{
        self.polygon == other.polygon
    }
}


pub struct GridNetwork{
    pub cells:Vec<GridCell>
}

impl GridNetwork{
    pub fn new(cells:Vec<GridCell>)->GridNetwork{
        GridNetwork{cells:cells}
    }
    pub fn query_node(self:&GridNetwork,node: &coords::Coordinate)->Vec<&GridCell>{
        // returns a copy of the cells that contain the node
        let mut cells:Vec<&GridCell> = Vec::new();
        for cell in self.cells.iter(){
            if cell.polygon.nodes.contains(&node){
                cells.push(cell);
            }
        }
        cells
    }
    pub fn query_edge(self:&GridNetwork,edge:&coords::Edge)->Vec<&GridCell>{
        // returns a copy of the cells that contain the edge
        let mut cells:Vec<&GridCell> = Vec::new();
        for cell in self.cells.iter(){
            if cell.polygon.to_edges().contains(&edge){
                cells.push(cell);
            }
        }
        cells
    }
    pub fn query_neighbors(self:&GridNetwork,cell: &GridCell)->Vec<&GridCell>{
        // returns references to the cells that are neighbors of the cell
        let mut cells:Vec<&GridCell> = Vec::new();
        for edge in cell.polygon.to_edges(){
            let neighbors = self.query_edge(&edge);
            for neighbor in neighbors{
                if neighbor != cell{
                    cells.push(neighbor);
                }
            }
        }
        cells
    }
    pub fn query_cell(self:&GridNetwork,cell:GridCell)->Result<&GridCell,()> {
        for _cell in self.cells.iter() {
            if cell == *_cell {
                return Result::Ok(_cell);
            }
        }
        Result::Err(())
    }
    pub fn max_value(self:&GridNetwork)->f64{
        let mut max:f64 = 0.0;
        for cell in self.cells.iter(){
            if cell.value > max{
                max = cell.value;
            }
        }
        max
    }
    pub fn min_value(self:&GridNetwork)->f64{
        let mut min:f64 = std::f64::MAX;
        for cell in self.cells.iter(){
            if cell.value < min{
                min = cell.value;
            }
        }
        min
    }
    pub fn average_value(self:&GridNetwork)->f64{
        let mut sum:f64 = 0.0;
        for cell in self.cells.iter(){
            sum += cell.value;
        }
        sum / self.cells.len() as f64
    }
}
