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
        let [a, b, c, d] = crate::basis::make_orthonormal_basis(axis_1, axis_2)?;

        Some(Self::from_orthonormal_basis(
            a, b, c, d, angle, 0.0,
        ))
    }
    pub fn from_rotation_arc(from: Vec4, to: Vec4) -> Self {
        Self::from_axes_angle(from, to, from.normalize().dot(to.normalize()).acos()).unwrap_or(Rot4::IDENTITY)
    }

    pub fn from_rotation_matrix(mat: Mat4) -> Self {
        let (left, right) = factor_cayley(mat);
        Self(left, right)
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

/// Decomposes a 4D rotation into two isoclinic rotations represented by
/// unit quaternions.
// https://digital.csic.es/bitstream/10261/132980/1/ON-CAYLEYS.pdf
fn factor_cayley(m: Mat4) -> (Quat, Quat) {
    let a00 = m.x_axis.x;
    let a10 = m.x_axis.y;
    let a20 = m.x_axis.z;
    let a30 = m.x_axis.w;
    let a01 = m.y_axis.x;
    let a11 = m.y_axis.y;
    let a21 = m.y_axis.z;
    let a31 = m.y_axis.w;
    let a02 = m.z_axis.x;
    let a12 = m.z_axis.y;
    let a22 = m.z_axis.z;
    let a32 = m.z_axis.w;
    let a03 = m.w_axis.x;
    let a13 = m.w_axis.y;
    let a23 = m.w_axis.z;
    let a33 = m.w_axis.w;

    let w = -a00 - a11 - a22 - a33;
    let s10 = -a10 + a01;
    let s12 = -a12 + a21;
    let s13 = -a31 + a13;
    let s20 = -a20 + a02;
    let s23 = -a23 + a32;
    let s30 = -a30 + a03;

    // `w` is the real component.
    let rr = Quat::from_xyzw(s23 + s10, s13 + s20, s12 + s30, w).normalize();
    let lr = Quat::from_xyzw(-s23 + s10, -s13 + s20, -s12 + s30, w).normalize();

    (lr, rr)
}
