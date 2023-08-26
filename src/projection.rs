use glam::{Mat4, Vec3, Vec4};

/// A `4D -> 3D` stereographic projection from `S3 -> R3`.
#[derive(Copy, Clone, Debug)]
pub struct Projection {
    /// Center of projection
    pub p: Vec4,
    /// Three elements of the orthonormal basis for the tangent plane
    /// of a point -- plane is centered at 0 rather than p.
    pub local_x: Vec4,
    pub local_y: Vec4,
    pub local_z: Vec4,
    /// Inverse matrix, used when projecting to strip a dimension out.
    pub matrix_inverse: Mat4,
}

impl Projection {
    /// The unit projection uses [`Vec4::W`] as the center of projection.
    pub const UNIT_PROJECTION: Projection = Projection {
        p: Vec4::W,
        local_x: Vec4::X,
        local_y: Vec4::Y,
        local_z: Vec4::Z,
        matrix_inverse: Mat4::IDENTITY,
    };

    /// Creates a new projection with `p` as the center.
    ///
    /// Picks an arbitrary orthonormal basis for `local_x`, `local_y`, and `local_z`.
    pub fn new(p: Vec4) -> Self {
        let q = crate::basis::any_orthogonal_vector_to_vector(p);
        let [p, q, r, s] = crate::basis::make_orthonormal_basis(p, q)
            .expect("Center of projection cannot be zero.");
        Self::from_orthonormal_basis(p, q, r, s)
    }

    /// Constructs a projection from an existing orthonormal basis.
    pub fn from_orthonormal_basis(p: Vec4, q: Vec4, r: Vec4, s: Vec4) -> Self {
        let mat = Mat4::from_cols(q, r, s, p).transpose();
        Self {
            p,
            local_x: q,
            local_y: r,
            local_z: s,
            matrix_inverse: mat,
        }
    }

    /// Projects a point according to the stereographic projection described
    /// by `self`.
    pub fn project(self, point: Vec4) -> Vec3 {
        let x = point;
        let sum = x + self.p;
        let projected = 2.0 * ((sum / self.p.dot(sum)) - self.p);
        // When multiplied by `matrix_inverse`, this will make the `w` component 0.
        (self.matrix_inverse * projected).truncate()
    }

    /// Projects a tangent vector to a point in an angle-preserving way.
    ///
    /// If:
    /// - `normal` is tangent to the hypersphere at `at`, ie `normal.dot(at) == 0`
    /// - `at` is part of some shape `K` which lives on the surface of the hypersphere
    /// - `K'` is the projection of `K`
    /// - `at'` is the projection of `at`
    /// - `x(t) = (at + t * normal).normalize()` is a spherical line segment for `t > 0`
    /// - `x(t)` intersects `K` at an angle `alpha`
    /// - `n'` is `self.project_normal(at, normal)`
    ///
    /// Then:
    /// - `y(t) = at' + t * n'` intersects `K'` at an angle `alpha`.
    ///
    /// This is meaningful because it means that normal vectors remain normal vectors.
    pub fn project_normal(self, at: Vec4, normal: Vec4) -> Vec3 {
        let xp = at.dot(self.p) + 1.0;
        let np = normal.dot(self.p);

        let projected = normal * xp - np * (at + self.p);
        (self.matrix_inverse * projected).truncate().normalize()
    }
}
