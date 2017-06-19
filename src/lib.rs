pub struct Point(f64, f64);
pub struct PointWithValue(Point, f64);
pub type Polygon = Vec<PointWithValue>;

enum GridValue {
    Exterior,
    Boundary(f64),
    Interior(f64),
}

struct Grid {
    level: u8,
    grid: Vec<GridValue>,
}

struct Transform {
    offset: Point,
    scaling: f64,
}

pub struct HarmonicMap {
    grid: Grid,
    transform: Transform,
}

pub fn init(polygon: Polygon, levels: u8) -> HarmonicMap {
    // TODO
    let grid = Grid { level: levels, grid: Vec::new() };
    let transform = Transform { offset: Point(0.0, 0.0), scaling: 0.0 };
    HarmonicMap { grid: grid, transform: transform }
}

pub fn eval(map: HarmonicMap, p: Point) -> f64 {
    // TODO
    0.0
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
    }
}
