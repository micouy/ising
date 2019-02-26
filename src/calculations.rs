//! Utilities for calculations and measurements.

/// Calculates average energy fluctuation at the given temperature from a vector
/// of energy levels.
pub fn calc_dE(Es: &[u64], T: f64) -> f64 {
    let n = Es.len() as f64;
    let avg_E_sq = (Es.iter().fold(0, |sum, E| sum + E.pow(2)) as f64) / n;
    let avg_E = (Es.iter().sum::<u64>() as f64) / n;

    (avg_E_sq - avg_E.powi(2)) / T
}

/// Calculates energy fluctuations from a vector of magnetization levels.
pub fn calc_X(Is: &[f64]) -> f64 {
    let n = Is.len() as f64;
    let avg_I_sq = (Is.iter().fold(0.0, |sum, I| sum + I.powi(2)) as f64) / n;
    let avg_I = (Is.iter().sum::<f64>() as f64) / n;

    avg_I_sq - avg_I.powi(2)
}

/// Calculates the probability of a flip based on the energy difference it would
/// cause and the temperature.
pub fn calc_flip_probability(E_diff: f64, K: f64, T: f64) -> f64 {
    if E_diff < 0.0 {
        1.0
    } else {
        std::f64::consts::E.powf(-E_diff / (K * T))
    }
}

#[cfg(test)]
mod test {
    use ::pretty_assertions::assert_eq;

    use super::*;

    fn float_error(x: f64, t: f64) -> f64 {
        (x - t).abs() / t
    }

    #[test]
    fn test_calculate_energy_fluctuation() {
        let Es = &[3, 5, 10, 2];
        let T = 0.5;

        let dE = calc_dE(Es, T);

        assert!(float_error(dE, 19.0) < 0.01);
    }

    #[test]
    fn test_caluculate_magnetic_susceptibility() {
        let Is = &[0.2, 0.4, 0.6, 0.8];

        let X = calc_X(Is);

        assert!(float_error(X, 0.05) < 0.01);
    }

    #[test]
    fn test_calculate_flip_probability() {
        let k = 1.0;
        let T = 10.0;

        let E_diff = -10.0;
        let probability = calc_flip_probability(E_diff, k, T);

        assert!(float_error(probability, 1.0) < 0.01);

        let E_diff = 10.0;
        let probability = calc_flip_probability(E_diff, k, T);

        assert!(float_error(probability, 0.37) < 0.01);
    }
}
