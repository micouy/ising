//! Utilities for calculations and measurements.

fn calc_fluctuations(Es: Vec<u64>) -> f64 {
    let n = Es.len() as f64;
    let avg_E_sq = (Es.iter().fold(0, |sum, E| sum + E.pow(2)) as f64) / n;
    let avg_E = (Es.iter().fold(0, |sum, E| sum + E) as f64) / n;

    avg_E_sq - avg_E.powi(2)
}

#[cfg(test)]
mod test {
    use ::pretty_assertions::assert_eq;
    use super::*;

    #[test]
    fn test_calculate_fluctuations() {
        let Es = vec![3, 5, 10, 2];

        let fluctuations = calc_fluctuations(Es);

        assert_eq!(fluctuations, 9.5);
    }
}
