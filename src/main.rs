// testing harness

use glam::{Mat4, Vec4};
use hypersphere::rotation::Rot4;

fn almost_eq(p: Vec4, q: Vec4) {
    if (p - q).length() > 0.0001 {
        panic!("Distance {:?} too large!", p - q);
    }
}

fn main() {
    let pts = (-7..=7).map(|x| x as f32)
        .flat_map(move |x| {
            (-7..=7).map(|x| x as f32)
                .flat_map(move |y| {
                    (-7..=7).map(|x| x as f32)
                        .flat_map(move |z| {
                            (-7..=7).map(|x| x as f32).map(move |w| Vec4::new(x, y, z, w).normalize())
                        })
                })
        })
        .filter(|x| !x.is_nan())
        .collect::<Vec<_>>();
    // println!("{:?}", pts);

    {
        let p = Vec4::new(-0.70710677, -0.70710677, 0.0, 0.0);
        let q = Vec4::new(-0.74740934, -0.66436386, 0.0, 0.0);

        let mut axes = [p, q, Vec4::X, Vec4::Y, Vec4::Z, Vec4::W];
        hypersphere::gram_schmidt::gram_schmidt(&mut axes);
        println!("[");
        for a in axes {
            println!("    {:?},", a);
        }
        println!("]");
        make_mat4_rot(p, q);
    }

    for i in 0..pts.len() {
        let p = pts[i];
        for j in (i + 1)..pts.len() {
            let q = pts[j];
            if (p + q).length_squared() < 0.00000001 {
                continue;
            }
            // println!("Checking {:?}, {:?}", p, q);
            let rot1 = Rot4::from_rotation_arc(p, q);
            // let mat = make_mat4_rot(p, q);
            almost_eq(rot1.mul_vec4(p), q);
            // almost_eq(mat.mul_vec4(p), q);
            almost_eq(rot1.inverse().mul_vec4(q), p);
            //
            //         // for &r in &pts {
            //         //     let rot2 = Rot4::from_rotation_arc(q, r);
            //         //     almost_eq(rot2.mul_vec4(rot1.mul_vec4(p)), r);
            //         //     almost_eq(rot2.mul_rot4(rot1).mul_vec4(p), r);
            //         // }
        }
    }

    // for &p in &pts {
    //     for &q in &pts {
    //         if (p + q).length_squared() < 0.00000001 {
    //             continue;
    //         }
    //         // println!("Checking {:?}, {:?}", p, q);
    //         // let rot1 = Rot4::from_rotation_arc(p, q);
    //         let mat = make_mat4_rot(p, q);
    //     }
    // }
}

fn make_mat4_rot(from: Vec4, to: Vec4) -> Mat4 {
    let angle = from.normalize().dot(to.normalize()).acos();

    let Some([a, b, c, d]) = hypersphere::gram_schmidt::orthonormal_basis(from, to) else {
        return Mat4::IDENTITY;
    };

    let orig_mat = Mat4::from_cols(a, b, c, d);
    {
        let prod = orig_mat * orig_mat.transpose();
        almost_eq(prod.x_axis, Vec4::X);
        almost_eq(prod.y_axis, Vec4::Y);
        almost_eq(prod.z_axis, Vec4::Z);
        almost_eq(prod.w_axis, Vec4::W);
    }
    let (sin_ab, cos_ab) = angle.sin_cos();
    let (sin_cd, cos_cd) = (0.0f32).sin_cos();
    let rot_mat = Mat4::from_cols(
        Vec4::new(cos_ab, sin_ab, 0.0, 0.0),
        Vec4::new(-sin_ab, cos_ab, 0.0, 0.0),
        Vec4::new(0.0, 0.0, cos_cd, sin_cd),
        Vec4::new(0.0, 0.0, -sin_cd, cos_cd),
    );
    let prod_mat = orig_mat * rot_mat * orig_mat.transpose();
    prod_mat
}
