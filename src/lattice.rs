//! Stuff related to the spin lattice.

#![allow(non_snake_case)]
use ::ndarray::prelude::*;
use ::rand::prelude::*;

/// Struct containing the spin lattice and its size.
pub struct Lattice {
    size: usize,
    inner: Array2<i8>,
}

impl Lattice {
    /// Create a new [`Lattice`] of a certain size with randomly generated
    /// spins.
    pub fn new(size: usize) -> Self {
        let mut rng = thread_rng();
        let spins: [i8; 2] = [-1, 1];
        let inner = Array2::from_shape_fn((size, size), |_| {
            *spins[..].choose(&mut rng).unwrap()
        });

        Self { size, inner }
    }

    /// Create a new `Lattice` from `Array2<i8>`.
    ///
    /// # Examples
    ///
    /// ```
    /// # fn main() -> Result<(), Box<std::error::Error>> {
    /// # use ::ndarray::prelude::*;
    /// # use ising_lib::Lattice;
    /// let array = Array::from_shape_vec((2, 2), vec![1, -1, 1, -1])?;
    /// let lattice = Lattice::from_array(array);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Panics
    ///
    /// The function **must** panic if `array` is not
    /// [`square`][ndarray::ArrayBase::is_square] or if any of the spins is neither `-1` nor `1`.
    ///
    /// ```should_panic
    /// # fn main() -> Result<(), Box<std::error::Error>> {
    /// # use ::ndarray::prelude::*;
    /// # use ising_lib::Lattice;
    /// //                                             ↓ incorrect spin value
    /// let array = Array::from_shape_vec((2, 2), vec![5, -1, 1, -1])?;
    /// let lattice = Lattice::from_array(array);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// ```should_panic
    /// # fn main() -> Result<(), Box<std::error::Error>> {
    /// # use ::ndarray::prelude::*;
    /// # use ising_lib::Lattice;
    /// //                                 ↓  ↓ array isn't square
    /// let array = Array::from_shape_vec((1, 4), vec![1, 1, 1, 1])?;
    /// let lattice = Lattice::from_array(array);
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_array(array: Array2<i8>) -> Self {
        assert!(array.is_square(), "Array is not square.");
        assert!(
            array.iter().all(|spin| *spin == 1 || *spin == -1),
            "Invalid spin value."
        );

        Lattice {
            size: array.shape()[0],
            inner: array,
        }
    }

    fn size(&self) -> usize {
        self.size
    }

    fn roll_index(&self, ix: usize, amt: isize) -> usize {
        let size = self.size as isize;

        ((ix as isize + size + amt) % size) as usize
    }

    fn sum_all_neighbors(&self, (i, j): (usize, usize)) -> i8 {
        assert!(i < self.size && j < self.size);

        [
            (self.roll_index(i, 1), j),
            (self.roll_index(i, -1), j),
            (i, self.roll_index(j, 1)),
            (i, self.roll_index(j, -1)),
        ]
        .iter()
        .map(|ix| self.inner[*ix])
        .sum()
    }

    fn sum_two_neighbors(&self, (i, j): (usize, usize)) -> i8 {
        assert!(i < self.size && j < self.size);

        [(self.roll_index(i, 1), j), (i, self.roll_index(j, 1))]
            .iter()
            .map(|ix| self.inner[*ix])
            .sum()
    }

    /// Calculate the energy difference the flip of the spin would cause.
    /// Used to determine the probability of a flip.
    ///
    /// Lattice before:          Lattice after:
    /// ##| a|##                 ##| a|##
    /// --------                 --------
    ///  b| S| c                  b|-S| c
    /// --------                 --------
    /// ##| d|##                 ##| d|##
    ///
    /// d_E = E_2 - E_1 =
    ///  = ((-J) * (-S) * (a + b + c + d)) - ((-J) * S * (a + b + c + d)) =
    ///  = -J * (a + b + c + d) * ((-S) - S) =
    ///  = -2 * -J * (a + b + c + d) * S =
    ///  = 2 * J * S * (a + b + c + d)
    fn calc_dE(&self, (i, j): (usize, usize), J: f32) -> f32 {
        assert!(i < self.size && j < self.size, "Index out of bounds.");

        2.0 * J * ((self.inner[(i, j)] * self.sum_all_neighbors((i, j))) as f32)
    }

    fn calc_E(&self, J: f32) -> f32 {
        self.inner
            .indexed_iter()
            .map(|(ix, s)| -J * (s * self.sum_two_neighbors(ix)) as f32)
            .sum()
    }
}

#[cfg(test)]
mod test {
    use ::pretty_assertions::assert_eq;

    use super::*;

    #[test]
    fn test_create_lattice() {
        let t_size = 40;
        let lattice = Lattice::new(t_size);

        assert_eq!(lattice.size(), t_size);
    }

    #[test]
    fn test_create_lattice_from_array() {
        let t_size = 2;
        let t_array =
            Array::from_shape_vec((t_size, t_size), vec![1, -1, 1, -1])
                .unwrap();
        let lattice = Lattice::from_array(t_array);

        assert_eq!(lattice.size(), t_size);
    }

    #[test]
    #[should_panic]
    fn test_panic_on_create_lattice_from_invalid_array() {
        let t_invalid_spin_array =
            Array::from_shape_vec((2, 2), vec![5, -1, 1, -1]).unwrap();
        let _ = Lattice::from_array(t_invalid_spin_array);

        let t_invalid_dimensions_array =
            Array::from_shape_vec((1, 4), vec![1, 1, 1, 1]).unwrap();
        let _ = Lattice::from_array(t_invalid_dimensions_array);
    }

    #[test]
    fn test_sum_neighbors() {
        let t_size = 3;
        let spins = [-1, -1, 1, 1, 1, 1, 1, 1, -1];
        let t_array =
            Array::from_shape_vec((t_size, t_size), spins.to_vec()).unwrap();
        let lattice = Lattice::from_array(t_array);

        let sum = lattice.sum_all_neighbors((1, 1));
        let t_sum = -1 + 1 + 1 + 1;

        assert_eq!(sum, t_sum);
    }

    #[test]
    fn test_calculate_cells_dE() {
        let t_array =
            Array::from_shape_vec((3, 3), vec![-1, -1, 1, 1, 1, 1, -1, 1, 1])
                .unwrap();
        let lattice = Lattice::from_array(t_array);
        let J = 1.0;
        let dE = lattice.calc_dE((1, 1), J);
        let t_dE =
            2.0 * J * 1 as f32 * lattice.sum_all_neighbors((1, 1)) as f32;

        assert_eq!(dE, t_dE);
    }

    #[test]
    fn test_caluclate_E() {
        let t_size = 2;
        let t_array =
            Array::from_shape_vec((t_size, t_size), vec![-1, -1, 1, 1])
                .unwrap();
        let lattice = Lattice::from_array(t_array.clone());
        let J = 1.0;

        let E = lattice.calc_E(J);
        let t_E = 0.0;

        assert_eq!(E, t_E);
    }
}
