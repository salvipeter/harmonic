#[derive(Clone,Copy,Debug,PartialEq)]
pub struct Point(f64, f64);

pub struct PointWithValue(Point, f64);

pub type Polygon = Vec<PointWithValue>;

struct Pixel(usize, usize);

#[derive(Clone)]
enum GridValue {
    Untyped,
    Exterior,
    Boundary(f64),
    Interior(f64),
}

struct Grid {
    level: usize,
    data: Vec<GridValue>,
}

struct Transform {
    offset: Point,
    scaling: f64,
}

pub struct HarmonicMap {
    grid: Grid,
    transform: Transform,
}

fn bounding_box(polygon: Polygon) -> (Point, Point) {
    let mut min = polygon.first().unwrap().0;
    let mut max = min;
    for p in polygon {
        let q = p.0;
        if q.0 < min.0 { min.0 = q.0 }
        if q.1 < min.1 { min.1 = q.1 }
        if q.0 > max.0 { max.0 = q.0 }
        if q.1 > max.1 { max.1 = q.1 }
    }
    (min, max)
}

fn fromPixel(transform: Transform, p: Pixel) -> Point {
    let Point(x0, y0) = transform.offset;
    Point(x0 + p.0 as f64 * transform.scaling,
          y0 + p.1 as f64 * transform.scaling)
}

fn toPixel(transform: Transform, p: Point) -> Pixel {
    let Point(x0, y0) = transform.offset;
    Pixel((x0 + p.0 * transform.scaling).round() as usize,
          (y0 + p.1 * transform.scaling).round() as usize)
}

pub fn init(polygon: Polygon, levels: usize) -> HarmonicMap {
    let n = 2^levels;

    let mut grid = Grid { level: levels, data: Vec::new() };
    grid.data.resize(n * n, GridValue::Untyped);

    let (min, max) = bounding_box(polygon);
    let len = (max.0 - min.0).max(max.1 - min.1);
    let transform = Transform { offset: min, scaling: len / n as f64 };

    HarmonicMap { grid: grid, transform: transform }
}

pub fn eval(map: HarmonicMap, p: Point) -> f64 {
    // TODO
    0.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::bounding_box;

    #[test]
    #[should_panic]
    fn bbox_empty() {
        bounding_box(Vec::new());
    }

    #[test]
    fn bbox() {
        assert_eq!(bounding_box(vec!(PointWithValue(Point(0.0, 0.0), 1.0),
                                     PointWithValue(Point(2.0, 5.0), 2.0),
                                     PointWithValue(Point(8.0, 3.0), 5.0),
                                     PointWithValue(Point(-1.0, -3.0), 6.0))),
                   (Point(-1.0,-3.0), Point(8.0,5.0)));
    }
}
