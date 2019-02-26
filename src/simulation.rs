//! Utitilities for simulation.

#[cfg(test)]
mod test {
    use ::pretty_assertions::assert_eq;

    use super::*;

    struct TParams {
        t: i32,
    }

    #[test]
    fn test_simulate() {
        let t_params = TParams { t: 0 };
        let t_setup = |params: TParams| -> i32 { params.t };
        let t_step = |i: usize, state: &mut i32| -> bool {
            *state += 1;

            *state < 10
        };

        let state = simulate(t_params, t_setup, t_step);

        assert_eq!(state, 10);
    }
}
