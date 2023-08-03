// Credit: https://www.qfbox.info/epermute
/// Computes the even permutations of `values`.
pub fn even_permutations<T: PartialOrd>(values: &mut [T], mut inspector: impl FnMut(&[T])) {
    let mut parity_flip = false;

    values.sort_by(|x, y| x.partial_cmp(y).unwrap());

    inspector(values);

    'outer: loop {
        if let Some(i) = (0..values.len() - 1)
            .filter(|&x| values[x] < values[x + 1])
            .max()
        {
            let j = (i..values.len())
                .filter(|&x| values[i] < values[x])
                .max()
                .unwrap();
            values.swap(i, j);
            parity_flip = !parity_flip;
            values[i + 1..].reverse();
            parity_flip ^= (values.len() - (i + 1)) % 4 > 1;
        } else {
            values.reverse();
            parity_flip ^= values.len() % 4 > 1;
        }

        if !parity_flip {
            inspector(values);
        }

        for i in 0..values.len() - 1 {
            if values[i] < values[i + 1] {
                continue 'outer;
            }
        }

        return;
    }
}
