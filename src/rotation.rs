use crate::projection::Projection;
use glam::{Mat4, Quat, Vec4};

#[derive(Copy, Clone, Debug, Default)]
pub struct Rot4(Quat, Quat);

impl Rot4 {
    pub const IDENTITY: Self = Rot4(Quat::IDENTITY, Quat::IDENTITY);
    pub const NAN: Self = Rot4(Quat::NAN, Quat::NAN);
    pub const fn from_pq(p: Quat, q: Quat) -> Self {
        Rot4(p, q)
    }
    pub fn from_slice(slice: &[f32]) -> Self {
        assert!(slice.len() >= 8);
        Rot4(Quat::from_slice(slice), Quat::from_slice(&slice[4..]))
    }
    pub fn write_to_slice(self, slice: &mut [f32]) {
        assert!(slice.len() >= 8);
        self.0.write_to_slice(slice);
        self.1.write_to_slice(&mut slice[4..]);
    }
    pub fn from_orthonormal_basis(
        a: Vec4,
        b: Vec4,
        c: Vec4,
        d: Vec4,
        angle_ab: f32,
        angle_cd: f32,
    ) -> Self {
        // if Mat4::from_cols(a, b, c, d).determinant() < 0.0 {
        //     std::mem::swap(&mut c, &mut d);
        // }

        let orig_mat = Mat4::from_cols(a, b, c, d);
        let (sin_ab, cos_ab) = angle_ab.sin_cos();
        let (sin_cd, cos_cd) = angle_cd.sin_cos();
        let rot_mat = Mat4::from_cols(
            Vec4::new(cos_ab, sin_ab, 0.0, 0.0),
            Vec4::new(-sin_ab, cos_ab, 0.0, 0.0),
            Vec4::new(0.0, 0.0, cos_cd, sin_cd),
            Vec4::new(0.0, 0.0, -sin_cd, cos_cd),
        );
        let prod_mat = orig_mat * rot_mat * orig_mat.transpose();
        let (left, right) = factor_cayley(prod_mat);
        Rot4(left.normalize(), right.normalize())
    }
    pub fn from_axes_angle(axis_1: Vec4, axis_2: Vec4, angle: f32) -> Option<Self> {
        let [a, b, c, d] = crate::gram_schmidt::orthonormal_basis(axis_1, axis_2)?;

        Some(Self::from_orthonormal_basis(
            a, b, c, d, angle, 0.0,
        ))
    }
    pub fn from_rotation_arc(from: Vec4, to: Vec4) -> Self {
        Self::from_axes_angle(from, to, from.normalize().dot(to.normalize()).acos()).unwrap_or(Rot4::IDENTITY)
    }
    pub fn inverse(self) -> Self {
        Self(self.0.inverse(), self.1.inverse())
    }
    pub fn normalize(self) -> Self {
        Self(self.0.normalize(), self.1.normalize())
    }
    pub fn is_finite(self) -> bool {
        self.0.is_finite() && self.1.is_finite()
    }
    pub fn is_nan(self) -> bool {
        self.0.is_nan() || self.1.is_nan()
    }
    pub fn is_normalized(self) -> bool {
        self.0.is_normalized() && self.1.is_normalized()
    }
    pub fn is_near_identity(self) -> bool {
        (self.0.is_near_identity() && self.1.is_near_identity())
            || ((-self.0).is_near_identity() && (-self.1).is_near_identity())
    }
    pub fn slerp(self, end: Self, s: f32) -> Self {
        Self(self.0.slerp(end.0, s), self.1.slerp(end.1, s))
    }
    pub fn mul_vec4(self, rhs: Vec4) -> Vec4 {
        let result = self.0 * Quat::from_xyzw(rhs.y, rhs.z, rhs.w, rhs.x) * self.1;
        Vec4::new(result.w, result.x, result.y, result.z)
    }
    pub fn mul_rot4(self, rhs: Self) -> Self {
        Self(rhs.0 * self.0, self.1 * rhs.1)
        // Self(rhs.0 * self.0, rhs.1 * self.1)
    }
    pub fn mul_proj(self, rhs: Projection) -> Projection {
        Projection::from_orthonormal_basis(
            self.mul_vec4(rhs.p),
            self.mul_vec4(rhs.local_x),
            self.mul_vec4(rhs.local_y),
            self.mul_vec4(rhs.local_z),
        )
    }


    pub fn from_rotation_xy(angle: f32) -> Self {
        Self::from_orthonormal_basis(Vec4::X, Vec4::Y, Vec4::Z, Vec4::W, angle, 0.0)
    }
    pub fn from_rotation_xz(angle: f32) -> Self {
        Self::from_orthonormal_basis(Vec4::X, Vec4::Z, Vec4::Y, Vec4::W, angle, 0.0)
    }
    pub fn from_rotation_xw(angle: f32) -> Self {
        Self::from_orthonormal_basis(Vec4::X, Vec4::W, Vec4::Z, Vec4::Y, angle, 0.0)
    }
    pub fn from_rotation_yz(angle: f32) -> Self {
        Self::from_orthonormal_basis(Vec4::Y, Vec4::Z, Vec4::X, Vec4::W, angle, 0.0)
    }
    pub fn from_rotation_yw(angle: f32) -> Self {
        Self::from_orthonormal_basis(Vec4::Y, Vec4::W, Vec4::X, Vec4::Z, angle, 0.0)
    }
    pub fn from_rotation_zw(angle: f32) -> Self {
        Self::from_orthonormal_basis(Vec4::Z, Vec4::W, Vec4::X, Vec4::Y, angle, 0.0)
    }
}

fn quat_exp(x: Quat) -> Quat {
    let a = x.w;
    let v = x.xyz();
    let len_v = v.length();
    let norm = v.normalize_or_zero();
    let (sin, cos) = len_v.sin_cos();
    Quat::from_xyzw(norm.x * sin, norm.y * sin, norm.z * sin, cos) * a.exp()
}

fn mat_new_rows(
    a: f32, b: f32, c: f32, d: f32,
    e: f32, f: f32, g: f32, h: f32,
    i: f32, j: f32, k: f32, l: f32,
    m: f32, n: f32, o: f32, p: f32,
) -> Mat4 {
    Mat4::from_cols(
        Vec4::new(a, e, i, m),
        Vec4::new(b, f, j, n),
        Vec4::new(c, g, k, o),
        Vec4::new(d, h, l, p),
    )
}

// fn factor_isoclinic(m: Mat4) -> (Quat, Quat) {
//     let a00 = m.x_axis.x;
//     let a10 = m.x_axis.y;
//     let a20 = m.x_axis.z;
//     let a30 = m.x_axis.w;
//     let a01 = m.y_axis.x;
//     let a11 = m.y_axis.y;
//     let a21 = m.y_axis.z;
//     let a31 = m.y_axis.w;
//     let a02 = m.z_axis.x;
//     let a12 = m.z_axis.y;
//     let a22 = m.z_axis.z;
//     let a32 = m.z_axis.w;
//     let a03 = m.w_axis.x;
//     let a13 = m.w_axis.y;
//     let a23 = m.w_axis.z;
//     let a33 = m.w_axis.w;
//
//     let m = 0.25 * mat_new_rows(
//         a00 + a11 + a22 + a33, a10 - a01 - a32 + a23, a20 + a31 - a02 - a13, a30 - a21 + a12 - a03,
//         a10 - a01 + a32 - a23, -a00 - a11 + a22 + a33, a30 - a21 - a12 + a03, -a20 - a31 - a02 - a13,
//         a20 - a31 - a02 + a13, -a30 - a21 - a12 - a03, -a00 + a11 - a22 + a33, a10 + a01 - a32 - a23,
//         a30 + a21 - a12 - a03, a20 - a31 + a02 - a13, -a10 - a01 - a32 - a23, -a00 + a11 + a22 - a33,
//     );
//
// }

const A1: Mat4 = Mat4::from_cols(Vec4::W, Vec4::Z, Vec4::NEG_Y, Vec4::NEG_X);
const A2: Mat4 = Mat4::from_cols(Vec4::NEG_Z, Vec4::W, Vec4::X, Vec4::NEG_Y);
const A3: Mat4 = Mat4::from_cols(Vec4::Y, Vec4::NEG_X, Vec4::W, Vec4::NEG_Z);

fn factor_cayley(m: Mat4) -> (Quat, Quat) {
    let ma1 = m * A1;
    let ma2 = m * A2;
    let ma3 = m * A3;

    let l0m = -m + A1 * ma1 + A2 * ma2 + A3 * ma3;
    // let l1m = -0.25 * (ma1 + A1 * m + A3 * ma2 - A2 * ma3);
    // let l2m = -0.25 * (ma2 + A2 * m + A1 * ma3 - A3 * ma1);
    // let l3m = -0.25 * (ma3 + A3 * m + A2 * ma1 - A1 * ma2);
    let rr = l0m * l0m.x_axis.length_recip();
    // let rr = l0m * l0m.determinant().sqrt().sqrt().recip();
    let rl = m * rr.transpose();
    (
        Quat::from_xyzw(rl.x_axis.y, rl.x_axis.z, rl.x_axis.w, rl.x_axis.x),
        Quat::from_xyzw(rr.x_axis.y, rr.x_axis.z, rr.x_axis.w, rr.x_axis.x),
    )
}
