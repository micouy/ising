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
}
