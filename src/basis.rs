use glam::Vec4;

/// Picks some (possibly non-unit-length) vector orthogonal to both `u` and `v`.
///
/// If `u` and `v` are not linearly independent, this will yield [`Vec4::ZERO`].
pub fn any_orthogonal_vector_to_plane(u: Vec4, v: Vec4) -> Vec4 {
    let [x, y, z, w] = u.to_array();
    let [p, q, r, s] = v.to_array();

    let rwsz = r * w - s * z;
    let syqw = s * y - q * w;
    let qzry = q * z - r * y;

    let cross_x = Vec4::new(0.0, rwsz, syqw, qzry);
    if cross_x == Vec4::ZERO {
        let pwsx = p * w - s * x;
        let rxpz = r * x - p * z;

        let cross_y = Vec4::new(-rwsz, 0.0, pwsx, rxpz);

        if cross_y == Vec4::ZERO {
            let pyqx = p * y - q * x;
            let cross_z = Vec4::new(-syqw, -pwsx, 0.0, pyqx);

            if cross_z == Vec4::ZERO {
                Vec4::new(-qzry, -rxpz, -pyqx, 0.0)
            } else {
                cross_z
            }
        } else {
            cross_y
        }
    } else {
        cross_x
    }
}

/// Yields some (possibly non-unit-length) vector orthogonal to `v`.
/// The result will only be [`Vec4::ZERO`] if `v` is [`Vec4::ZERO`]
/// (modulo floating point error).
pub fn any_orthogonal_vector_to_vector(v: Vec4) -> Vec4 {
    let [x, y, z, w] = v.to_array();

    Vec4::new(w - z + y, z - w - x, w + x - y, y - z - x)
}

/// The 4D cross product.
///
/// Takes three vectors and yields a new vector that is linearly
/// independent to all three.
///
/// If the three vectors form a linearly dependent set, this will yield
/// [`Vec4::ZERO`].
pub fn cross_4d(u: Vec4, v: Vec4, w: Vec4) -> Vec4 {
    Vec4::new(
        -u.w * v.z * w.y + u.z * v.w * w.y + u.w * v.y * w.z - u.y * v.w * w.z - u.z * v.y * w.w + u.y * v.z * w.w,
        u.w * v.z * w.x - u.z * v.w * w.x - u.w * v.x * w.z + u.x * v.w * w.z + u.z * v.x * w.w - u.x * v.z * w.w,
        -u.w * v.y * w.x + u.y * v.w * w.x + u.w * v.x * w.y - u.x * v.w * w.y - u.y * v.x * w.w + u.x * v.y * w.w,
        u.z * v.y * w.x - u.y * v.z * w.x - u.z * v.x * w.y + u.x * v.z * w.y + u.y * v.x * w.z - u.x * v.y * w.z,
    )
}

/// Constructs an orthonormal basis for R4 from the plane spanned by `u` and `v`.
///
/// If `p` and `q` are not linearly independent, this will return `None`.
///
/// When this returns `Some([a, b, c, d])`, `a` is guaranteed to be `p.normalize()`,
/// however `b` could point in the "opposite direction" than `q` in the plane of `p q`.
pub fn make_orthonormal_basis(p: Vec4, q: Vec4) -> Option<[Vec4; 4]> {
    let r = any_orthogonal_vector_to_plane(p, q);
    let s = cross_4d(p, q, r);
    let q = cross_4d(p, r, s);
    Some([
        p.try_normalize()?,
        q.try_normalize()?,
        r.try_normalize()?,
        s.try_normalize()?,
    ])
}
