#![allow(non_snake_case)]

use ::ndarray::prelude::*;
use ::rand::prelude::*;

#[derive(Debug)]
enum LatticeCreationError {
    UnequalDimensions,
    InvalidSpinValue,
}

struct Lattice {
    size: usize,
    inner: Array2<i8>,
}

impl Lattice {
    fn new(size: usize) -> Self {
        let mut rng = thread_rng();
        let spins: [i8; 2] = [-1, 1];
        let inner = Array2::from_shape_fn((size, size), |_| {
            *spins[..].choose(&mut rng).unwrap()
        });

        Self { size, inner }
    }

    fn try_from_array(array: Array2<i8>) -> Result<Self, LatticeCreationError> {
        if array.shape()[0] != array.shape()[1] {
            return Err(LatticeCreationError::UnequalDimensions);
        }

        if !array.iter().all(|spin| spin.abs() == 1) {
            return Err(LatticeCreationError::InvalidSpinValue);
        }

        Ok(Lattice { size: array.shape()[0], inner: array })
    }

    fn size(&self) -> usize {
        self.size
    }

    /// Calculate the energy difference the flip of the spin would cause.
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

        let neighbours = [
            ((i + 1) % self.size, j),
            ((i - 1) % self.size, j),
            (i, (j + 1) % self.size),
            (i, (j + 1) % self.size),
        ];

        let neighbours_sum = neighbours
            .iter()
            .map(|ix| self.inner[*ix])
            .fold(0, |sum, n| sum + n);

        let s = self.inner[(i, j)];

        2.0 * J * ((s * neighbours_sum) as f32)
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
        let t_array = Array::from_shape_vec(
            (t_size, t_size),
            vec![1, -1, 1, -1]
        ).unwrap();
        let lattice = Lattice::try_from_array(t_array).unwrap();

        assert_eq!(lattice.size(), t_size);
    }

    #[test]
    fn test_error_on_create_lattice_from_invalid_array() {
        let t_invalid_spin_array = Array::from_shape_vec(
            (2, 2),
            vec![5, -1, 1, -1]
        ).unwrap();
        let result = Lattice::try_from_array(t_invalid_spin_array);

        assert!(result.is_err());

        let t_invalid_dimensions_array = Array::from_shape_vec(
            (1, 4),
            vec![1, 1, 1, 1]
        ).unwrap();
        let result = Lattice::try_from_array(t_invalid_dimensions_array);

        assert!(result.is_err());
    }

    #[test]
    fn test_calculate_cells_dE() {
        let t_array = Array::from_shape_vec(
            (3, 3),
            vec![
                -1, -1,  1,
                 1,  1,  1,
                -1,  1,  1,
            ]
        ).unwrap();
        let lattice = Lattice::try_from_array(t_array).unwrap();
        let J = 1.0;
        let dE = lattice.calc_dE((1, 1), J);
        let t_dE = 2.0 * J * 1 as f32 * (-1 + 1 + 1 + 1) as f32;

        assert_eq!(dE, t_dE);
    }
}
