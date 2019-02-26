//! Utilities for calculations and measurements.

fn calc_fluctuations(Es: Vec<u64>, T: f64) -> f64 {
    let n = Es.len() as f64;
    let avg_E_sq = (Es.iter().fold(0, |sum, E| sum + E.pow(2)) as f64) / n;
    let avg_E = (Es.iter().sum::<u64>() as f64) / n;

    (avg_E_sq - avg_E.powi(2)) / T
}

fn calc_mag_susceptibility(Is: Vec<f64>) -> f64 {
    let n = Is.len() as f64;
    let avg_I_sq = (Is.iter().fold(0.0, |sum, I| sum + I.powi(2)) as f64) / n;
    let avg_I = (Is.iter().sum::<f64>() as f64) / n;

    avg_I_sq - avg_I.powi(2)
}

#[cfg(test)]
mod test {
    use ::pretty_assertions::assert_eq;

    use super::*;

    fn float_error(x: f64, t: f64) -> f64 {
        (x - t).abs() / t
    }

    #[test]
    fn test_calculate_fluctuations() {
        let Es = vec![3, 5, 10, 2];
        let T = 0.5;

        let fluctuations = calc_fluctuations(Es, T);

        assert!(float_error(fluctuations, 19.0) < 0.01);
    }

    #[test]
    fn test_caluculate_magnetic_susceptibility() {
        let Is = vec![0.2, 0.4, 0.6, 0.8];

        let mag_susceptibility = calc_mag_susceptibility(Is);

        assert!(float_error(mag_susceptibility, 0.05) < 0.01);
    }
}
