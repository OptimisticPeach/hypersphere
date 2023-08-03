use glam::Vec4;

pub fn gram_schmidt(vectors: &mut [Vec4]) {
    if vectors.len() == 0 {
        return;
    }

    vectors[0] = vectors[0].try_normalize().unwrap();
    for i in 1..vectors.len() {
        let projected: Vec4 = vectors[..i].iter().map(|x| x.dot(vectors[i]) * *x).sum();
        let diff = vectors[i] - projected;
        let length_sq = diff.length_squared();
        if length_sq < 0.0000000001 {
            vectors[i] = Vec4::ZERO;
        } else {
            vectors[i] = diff / length_sq.sqrt();
        }
    }

    return;
}
