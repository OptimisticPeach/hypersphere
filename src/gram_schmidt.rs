use glam::Vec4;

pub fn gram_schmidt(vectors: &mut [Vec4]) {
    if vectors.len() == 0 {
        return;
    }

    vectors[0] = vectors[0].try_normalize().unwrap();
    for i in 1..vectors.len() {
        let mut projected = vectors[i];
        for &v in &vectors[..i] {
            projected -= project(v, vectors[i]);
        }
        let length = projected.length();
        if length < 0.0001 {
            vectors[i] = Vec4::ZERO;
        } else {
            vectors[i] = projected.normalize_or_zero();
        }
        // println!("Projected: {projected:?}, len: {:?}", projected.length());
        // let projected: Vec4 = vectors[..i].iter().map(|x| x.dot(vectors[i]) * x.dot(*x).recip() * *x).sum();
        // let diff = vectors[i] - projected;
        // let length_sq = diff.length_squared();
        // if length_sq < 0.0000000001 {
        //     vectors[i] = Vec4::ZERO;
        // } else {
        //     vectors[i] = diff / length_sq.sqrt();
        // }
    }
}

pub fn orthonormal_basis(u: Vec4, v: Vec4) -> Option<[Vec4; 4]> {
    let mut axes = [u, v, Vec4::X, Vec4::Y, Vec4::Z, Vec4::W];
    gram_schmidt(&mut axes);
    if axes[1] == Vec4::ZERO {
        return None;
    }
    // println!("Axes Original: {:?}", axes);

    for i in 2..axes.len() {
        // let length = axes[i].length_squared();
        // println!("Checking i = {i}, len: {length}");
        if axes[i] != Vec4::ZERO {
            axes[2] = axes[i];
            // println!("i = {}", i);

            for j in i + 1..axes.len() {
                // let length = axes[j].length_squared();
                // println!("Checking j = {j}, len: {length}");
                if axes[j] != Vec4::ZERO {
                    axes[3] = axes[j];
                    // println!("i = {}, j = {}", i, j);
                    break;
                }
            }
            break;
        }
    }

    // let [a, b, mut c, mut d, _, _] = axes;
    // // println!("Axes: {} {} {} {}", a, b, c, d);
    // if Mat4::from_cols(a, b, c, d).determinant() < 0.0 {
    //     std::mem::swap(&mut c, &mut d);
    // }

    Some([axes[0], axes[1], axes[2], axes[3]])
}

fn project(u: Vec4, v: Vec4) -> Vec4 {
    let uurecip = u.dot(u);
    if uurecip == 0.0 {
        return Vec4::ZERO;
    }
    v.dot(u) * uurecip.recip() * u
}
