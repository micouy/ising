//! Stuff related to spin lattice.

#![allow(non_snake_case)]
use ::ndarray::prelude::*;
use ::rand::prelude::*;

/// Struct encapsulating the spin lattice and operations on it.
///
/// The lattice behaves like a torus - spins on the opposite edges are
/// considered each other's neighbors.
pub struct Lattice {
    size: usize,
    inner: Array2<i8>,
}

impl Lattice {
    /// Creates a new [`Lattice`] of a certain size with randomly generated
    /// spins.
    pub fn new(size: usize) -> Self {
        let mut rng = thread_rng();
        let spins: [i8; 2] = [-1, 1];
        let inner = Array2::from_shape_fn((size, size), |_| {
            *spins[..].choose(&mut rng).unwrap()
        });

        Self { size, inner }
    }

    /// Creates a new [`Lattice`] from [`Array2<i8>`][ndarray::Array2].
    ///
    /// # Examples
    ///
    /// ```
    /// # fn main() -> Result<(), Box<std::error::Error>> {
    /// # use ::ndarray::prelude::*;
    /// # use ising_lib::prelude::*;
    /// let array = Array::from_shape_vec((2, 2), vec![1, -1, 1, -1])?;
    /// let lattice = Lattice::from_array(array);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Panics
    ///
    /// The function will panic if `array` is not
    /// [`square`][ndarray::ArrayBase::is_square] or if any of the spins
    /// has incorrect value (neither `-1` nor `1`).
    ///
    /// ```should_panic
    /// # fn main() -> Result<(), Box<std::error::Error>> {
    /// # use ::ndarray::prelude::*;
    /// # use ising_lib::prelude::*;
    /// let array = Array::from_shape_vec((2, 2), vec![5, -1, 1, -1])?;
    /// //                                             ↑ incorrect spin value
    /// let lattice = Lattice::from_array(array);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// ```should_panic
    /// # fn main() -> Result<(), Box<std::error::Error>> {
    /// # use ::ndarray::prelude::*;
    /// # use ising_lib::prelude::*;
    /// let array = Array::from_shape_vec((1, 4), vec![1, 1, 1, 1])?;
    /// //                                 ↑  ↑ array isn't square
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

    /// Returns the size of the lattice.
    pub fn size(&self) -> usize {
        self.size
    }

    /// Returns the `(ith, jth)` spin value.
    fn get(&self, ix: (usize, usize)) -> i8 {
        self.inner[ix]
    }

    /// Returns the modulo operation on the `index + amt`, so that it stays
    /// within the bounds.
    fn roll_index(&self, (i, j): (usize, usize), amt: isize) -> (usize, usize) {
        let size = self.size as isize;

        (
            ((i as isize + size + amt) % size) as usize,
            ((j as isize + size + amt) % size) as usize,
        )
    }

    /// Returns the product of the `(ith, jth)` spin and the sum of all of its
    /// neighbors.
    fn spin_times_all_neighbors(&self, (i, j): (usize, usize)) -> i8 {
        assert!(i < self.size && j < self.size);

        let (i_1, j_1) = self.roll_index((i, j), -1);
        let (i_2, j_2) = self.roll_index((i, j), 1);

        [(i_1, j), (i_2, j), (i, j_1), (i, j_2)]
            .iter()
            .map(|ix| self.inner[*ix])
            .sum()
    }

    /// Returns the product of the `(ith, jth)` spin and the sum of two of its
    /// neighbors.
    fn spin_times_two_neighbors(&self, (i, j): (usize, usize)) -> i8 {
        assert!(i < self.size && j < self.size);

        let (i_r, j_r) = self.roll_index((i, j), 1);

        [(i_r, j), (i, j_r)].iter().map(|ix| self.inner[*ix]).sum()
    }

    /// Calculates the difference of energy that would be caused by
    /// flipping the `(ith, jth)` spin without actually doing it.
    /// Used to determine the probability of a flip.
    ///
    /// ```text
    /// Lattice before flip:     Lattice after flip:
    /// ##| a|##                 ##| a|##
    /// --------                 --------
    ///  b| s| c                  b|-s| c
    /// --------                 --------
    /// ##| d|##                 ##| d|##
    ///
    /// d_E = E_2 - E_1 =
    ///  = ((-J) * (-s) * (a + b + c + d)) - ((-J) * s * (a + b + c + d)) =
    ///  = -J * (a + b + c + d) * ((-s) - s) =
    ///  = -2 * -J * (a + b + c + d) * s =
    ///  = 2 * J * s * (a + b + c + d)
    /// ```
    ///
    /// # Panics
    ///
    /// The function will panic if the index is out of bounds.
    ///
    /// ```should_panic
    /// # use ising_lib::prelude::*;
    /// let lattice = Lattice::new(10);
    /// let _ = lattice.calc_dE((42, 0), 1.0);
    /// ```
    pub fn calc_dE(&self, (i, j): (usize, usize), J: f32) -> f32 {
        assert!(i < self.size && j < self.size);

        2.0 * J * f32::from(self.spin_times_all_neighbors((i, j)))
    }

    /// Calculates the energy of the lattice.
    pub fn calc_E(&self, J: f32) -> f32 {
        self.inner
            .indexed_iter()
            .map(|(ix, _)| -J * f32::from(self.spin_times_two_neighbors(ix)))
            .sum()
    }

    /// Calculates the magnetization of the lattice. The magnetization is
    /// a value in range `[0.0, 1.0]` and it is the absolute value of the mean
    /// spin value.
    pub fn calc_I(&self) -> f32 {
        f32::from(self.inner.sum().abs()) / self.size.pow(2) as f32
    }

    /// Flips the `(ith, jth)` spin.
    pub fn flip_spin(&mut self, (i, j): (usize, usize)) {
        assert!(i < self.size && j < self.size);

        *self.inner.get_mut((i, j)).unwrap() *= -1;
    }
}

#[cfg(test)]
mod test {
    use ::pretty_assertions::assert_eq;

    use super::*;

    #[test]
    fn test_create_lattice() {
        let lattice = Lattice::new(40);

        assert_eq!(lattice.size(), 40);
    }

    #[test]
    fn test_create_lattice_from_array() {
        let array = Array::from_shape_vec((2, 2), vec![1, -1, 1, -1]).unwrap();
        let lattice = Lattice::from_array(array);

        assert_eq!(lattice.size(), 2);
    }

    #[test]
    fn test_spin_times_neighbors() {
        let spins = [-1, -1, 1, 1, 1, 1, 1, 1, -1];
        let array = Array::from_shape_vec((3, 3), spins.to_vec()).unwrap();
        let lattice = Lattice::from_array(array);

        let product = lattice.spin_times_all_neighbors((1, 1));

        assert_eq!(product, 2);
    }

    #[test]
    fn test_calculate_dE() {
        let array =
            Array::from_shape_vec((3, 3), vec![-1, -1, 1, 1, 1, 1, -1, 1, 1])
                .unwrap();
        let lattice = Lattice::from_array(array);
        let J = 1.0;

        let dE = lattice.calc_dE((1, 1), J);

        assert_eq!(dE, 4.0);
    }

    #[test]
    fn test_caluclate_E() {
        let array = Array::from_shape_vec((2, 2), vec![-1, -1, 1, 1]).unwrap();
        let lattice = Lattice::from_array(array);
        let J = 1.0;

        let E = lattice.calc_E(J);

        assert_eq!(E, 0.0);
    }

    #[test]
    fn test_calculate_I() {
        let array = Array::from_shape_vec((2, 2), vec![-1, -1, -1, 1]).unwrap();
        let lattice = Lattice::from_array(array);

        let I = lattice.calc_I();

        assert_eq!(I, 0.5);
    }

    #[test]
    fn test_flip_spin() {
        let array = Array::from_shape_vec((2, 2), vec![-1, -1, -1, 1]).unwrap();
        let mut lattice = Lattice::from_array(array);

        lattice.flip_spin((1, 1));

        assert_eq!(lattice.get((1, 1)), -1);
    }
}
