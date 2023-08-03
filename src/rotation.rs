use crate::projection::Projection;
use glam::{Quat, Vec4};

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
        let [a, b, c, d] = [a, b, c, d].map(|x| Quat::from_vec4(Vec4::new(x.y, x.z, x.w, x.x)));

        let left_inner = b * a.conjugate() * (angle_ab + angle_cd) / 2.0;
        let right_inner = c.conjugate() * d * (angle_ab - angle_cd) / 2.0;

        let left = quat_exp(left_inner).normalize();
        let right = quat_exp(right_inner).normalize();

        Rot4(left, right)
    }
    pub fn from_axes_angle(axis_1: Vec4, axis_2: Vec4, angle: f32) -> Option<Self> {
        let mut axes = [axis_1, axis_2, Vec4::X, Vec4::Y, Vec4::Z, Vec4::W];
        super::gram_schmidt::gram_schmidt(&mut axes);
        if axes[1] == Vec4::ZERO {
            return None;
        }

        for i in 2..axes.len() - 2 {
            if axes[i] != Vec4::ZERO {
                axes[2] = axes[i];
                for j in i + 1..axes.len() - (i + 1) {
                    if axes[j] != Vec4::ZERO {
                        axes[3] = axes[j];
                        break;
                    }
                }
                break;
            }
        }

        Some(Self::from_orthonormal_basis(
            axes[0], axes[1], axes[2], axes[3], angle, 0.0,
        ))
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
    pub fn from_rotation_arc(from: Vec4, to: Vec4) -> Self {
        Self::from_axes_angle(from, to, from.dot(to).acos()).unwrap_or(Rot4::IDENTITY)
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
            self.mul_vec4(rhs.q),
            self.mul_vec4(rhs.r),
            self.mul_vec4(rhs.s),
        )
    }
}

fn quat_exp(x: Quat) -> Quat {
    let a = x.w;
    let v = x.xyz();
    let len_v = v.length();
    let norm = v / len_v;
    let (sin, cos) = len_v.sin_cos();
    Quat::from_xyzw(norm.x * sin, norm.y * sin, norm.z * sin, cos) * a.exp()
}
