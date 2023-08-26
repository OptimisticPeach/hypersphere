use glam::Vec4;

/// Creates the vertices of the 600 Cell.
///
/// Useful as a set of interesting debugging vertices.
///
/// This shape is also known sometimes as the hexachosichoron.
pub fn make_600_cell() -> Vec<Vec4> {
    let mut points = Vec::new();

    points.extend_from_slice(&[
        Vec4::X,
        Vec4::NEG_X,
        Vec4::Y,
        Vec4::NEG_Y,
        Vec4::Z,
        Vec4::NEG_Z,
        Vec4::W,
        Vec4::NEG_W,
    ]);

    for i in 0..16 {
        let s0 = i & 1 == 0;
        let s1 = i & 2 == 0;
        let s2 = i & 4 == 0;
        let s3 = i & 8 == 0;

        let sign = |bool_val| if bool_val { 0.5 } else { -0.5 };

        points.push(Vec4::new(sign(s0), sign(s1), sign(s2), sign(s3)));
    }

    let mut x = [
        (1.0 + 5.0f32.sqrt()) / 4.0,
        0.5,
        (1.0 + 5.0f32.sqrt()).recip(),
        0.0,
    ];

    crate::even_permutations::even_permutations(&mut x, |x| {
        let mut x_copied = [x[0], x[1], x[2], x[3]];
        let zero_idx: usize = if x[0] == 0.0 {
            0
        } else if x[1] == 0.0 {
            1
        } else if x[2] == 0.0 {
            2
        } else {
            3
        };

        x_copied.swap(zero_idx, 3);

        for parity in 0..8 {
            let s0 = parity & 1 == 0;
            let s1 = parity & 2 == 0;
            let s2 = parity & 4 == 0;

            let sign = |bool_val, val| -> f32 {
                if bool_val {
                    val
                } else {
                    -val
                }
            };

            let mut new_values = [
                sign(s0, x_copied[0]),
                sign(s1, x_copied[1]),
                sign(s2, x_copied[2]),
                0.0,
            ];
            new_values.swap(3, zero_idx);
            points.push(Vec4::from_slice(&new_values));
        }
    });

    points
}
