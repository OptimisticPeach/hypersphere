use glam::{Mat4, Vec3, Vec4};

#[derive(Copy, Clone, Debug)]
pub struct Projection {
    // Center of projection
    pub p: Vec4,
    // Three elements of the orthonormal basis for the tangent plane
    // of a point -- plane is centered at 0 rather than p.
    pub local_x: Vec4,
    pub local_y: Vec4,
    pub local_z: Vec4,
    pub matrix_inverse: Mat4,
}

impl Projection {
    pub const UNIT_PROJECTION: Projection = Projection {
        p: Vec4::W,
        local_x: Vec4::X,
        local_y: Vec4::Y,
        local_z: Vec4::Z,
        matrix_inverse: Mat4::IDENTITY,
    };

    pub fn new(p: Vec4) -> Self {
        let p = p
            .try_normalize()
            .expect("Center of projection cannot be zero.");
        let mut axes = [p, Vec4::X, Vec4::Y, Vec4::Z, Vec4::W];
        super::gram_schmidt::gram_schmidt(&mut axes);

        for i in 1..axes.len() - 1 {
            if axes[i] != Vec4::ZERO {
                axes[1] = axes[i];
                for j in i + 1..axes.len() - (i + 1) {
                    if axes[j] != Vec4::ZERO {
                        axes[2] = axes[j];
                        for k in j + 1..axes.len() - (j + 1) {
                            if axes[k] != Vec4::ZERO {
                                axes[3] = axes[k];
                                break;
                            }
                        }
                        break;
                    }
                }
                break;
            }
        }
        let [p, q, r, s, ..] = axes;
        Self::from_orthonormal_basis(p, q, r, s)
    }

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

    pub fn project(self, point: Vec4) -> Vec3 {
        let x = point;
        let sum = x + self.p;
        let projected = 2.0 * ((sum / self.p.dot(sum)) - self.p);
        (self.matrix_inverse * projected).truncate()
    }

    pub fn project_normal(self, at: Vec4, normal: Vec4) -> Vec3 {
        let xp = at.dot(self.p) + 1.0;
        let np = normal.dot(self.p);

        let projected = normal * xp - np * (at + self.p);
        (self.matrix_inverse * projected).truncate().normalize()
    }

    pub fn unproject_direction(self, direction: Vec3) -> Vec4 {
        direction.x * self.local_x + direction.y * self.local_y + direction.z * self.local_z
    }

    pub fn place_in_tangent_space(self, value: Vec3) -> Vec4 {
        (value.x * self.local_x + value.y * self.local_y + value.z * self.local_z + self.p).normalize()
    }
}
