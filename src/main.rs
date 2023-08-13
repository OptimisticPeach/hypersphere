// testing harness

use glam::Vec4;
use hexasphere::shapes::IcoSphere;
use hypersphere::projection::Projection;
use rand::{thread_rng, Rng};

fn main() {
    let ico = IcoSphere::new(30, |_| {});
    let points = ico.raw_points();

    let projection = Projection::new(Vec4::new(3.0, 5.0, -1.0, 7.0));
    let mut rng = thread_rng();

    let radii = points
        .iter()
        .map(|_| rng.gen_range(0.5f32..1.3))
        .collect::<Vec<_>>();

    let points = hypersphere::points_at(
        std::f32::consts::PI / 180.0,
        Vec4::new(1.0, 3.0, 4.0, 5.0),
        points,
        &radii,
    );

    // let mut points = hypersphere::sphere_at(
    //     5.0 * std::f32::consts::PI / 180.0,
    //     Vec4::new(1.0, 3.0, 4.0, 5.0).normalize(),
    //     points,
    // );

    println!(
        "Min: {:?}",
        points
            .iter()
            .map(|x| x.length())
            .min_by(|x, y| x.partial_cmp(y).unwrap())
    );

    println!(
        "Max: {:?}",
        points
            .iter()
            .map(|x| x.length())
            .max_by(|x, y| x.partial_cmp(y).unwrap())
    );

    // let projection = Projection::new(Vec4::new(3.0, 5.0, -1.0, 7.0));
    //
    // points.iter_mut()
    //     .for_each(|x| *x = projection.project(*x));
    //
    // println!(
    //     "Max w: {:?}",
    //     points
    //         .iter()
    //         .map(|x| x.w)
    //         .max_by(|x, y| x.partial_cmp(y).unwrap())
    // );
    //
    // println!(
    //     "Min w: {:?}",
    //     points
    //         .iter()
    //         .map(|x| x.w)
    //         .min_by(|x, y| x.partial_cmp(y).unwrap())
    // );
}
