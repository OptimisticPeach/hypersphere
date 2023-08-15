use crate::rotation::Rot4;
use glam::{Vec3, Vec3A, Vec4};

pub mod camera;
mod even_permutations;
mod gram_schmidt;
pub mod projection;
pub mod rotation;

pub fn make_points() -> Vec<Vec4> {
    let mut points = Vec::new();

    points.extend_from_slice(&[
        Vec4::X,
        Vec4::NEG_X,
        Vec4::Y,
        Vec4::NEG_Y,
        Vec4::Z,
        Vec4::NEG_Z,
        Vec4::W,
        Vec4::NEG_W,
    ]);

    for i in 0..16 {
        let s0 = i & 1 == 0;
        let s1 = i & 2 == 0;
        let s2 = i & 4 == 0;
        let s3 = i & 8 == 0;

        let sign = |bool_val| if bool_val { 0.5 } else { -0.5 };

        points.push(Vec4::new(sign(s0), sign(s1), sign(s2), sign(s3)));
    }

    let mut x = [
        (1.0 + 5.0f32.sqrt()) / 4.0,
        0.5,
        (1.0 + 5.0f32.sqrt()).recip(),
        0.0,
    ];

    even_permutations::even_permutations(&mut x, |x| {
        let mut x_copied = [x[0], x[1], x[2], x[3]];
        let zero_idx: usize = if x[0] == 0.0 {
            0
        } else if x[1] == 0.0 {
            1
        } else if x[2] == 0.0 {
            2
        } else {
            3
        };

        x_copied.swap(zero_idx, 3);

        for parity in 0..8 {
            let s0 = parity & 1 == 0;
            let s1 = parity & 2 == 0;
            let s2 = parity & 4 == 0;

            let sign = |bool_val, val| -> f32 {
                if bool_val {
                    val
                } else {
                    -val
                }
            };

            let mut new_values = [
                sign(s0, x_copied[0]),
                sign(s1, x_copied[1]),
                sign(s2, x_copied[2]),
                0.0,
            ];
            new_values.swap(3, zero_idx);
            points.push(Vec4::from_slice(&new_values));
        }
    });

    points
}

pub fn spherical_radius_to_size(angle: f32) -> f32 {
    angle.sin()
}

pub fn size_to_spherical_radius(size: f32) -> f32 {
    size.asin()
}

pub fn sphere_at(spherical_radius: f32, at: Vec4, points: &[Vec3A]) -> Vec<Vec4> {
    let at = at.normalize();
    let center = at * spherical_radius.cos();
    let size = spherical_radius_to_size(spherical_radius);
    let rotation = Rot4::from_rotation_arc(Vec4::W, at);
    points
        .iter()
        .map(|&x| rotation.mul_vec4((x * size).extend(0.0)) + center)
        .collect()
}

pub fn mesh_at(
    spherical_radius: f32,
    at: Vec4,
    points: &[Vec3],
    normals: &[Vec3],
    radii: &[f32],
) -> (Vec<Vec4>, Vec<Vec4>) {
    let at = at.normalize();
    let rotation = Rot4::from_rotation_arc(Vec4::W, at);
    let points = points
        .iter()
        .copied()
        .zip(radii.iter().copied())
        .map(|(pos, len)| {
            let radius = spherical_radius * len;
            let (size, center_scl) = radius.sin_cos();
            let center = at * center_scl;
            rotation.mul_vec4((pos * size).extend(0.0)) + center
        })
        .collect::<Vec<_>>();

    let normals = normals
        .iter()
        .copied()
        .zip(points.iter().copied())
        .map(|(normal, pos_4d)| {
            let rotated = rotation.mul_vec4(normal.extend(0.0));
            let gram_schmidt = rotated - rotated.dot(pos_4d) * pos_4d;
            gram_schmidt.normalize()
        })
        .collect::<Vec<_>>();

    (points, normals)
}
