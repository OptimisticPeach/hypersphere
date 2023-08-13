use buttery::Smoothed;
use crate::rotation::Rot4;

/// Represents 4D rotation interpolation through [`slerp`](Rot4::slerp).
pub struct Rotate4D;

impl Smoothed for Rotate4D {
    type Attribute = Rot4;
    fn drive(target: Rot4, current: Rot4, percent: f32) -> Rot4 {
        current.slerp(target, percent).normalize()
    }
}
