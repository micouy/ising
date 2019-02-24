//! Everything you need to run a simplified
//! [Ising model](https://en.wikipedia.org/wiki/Ising_model) simulation.
//! In this model every spin's value must be either `-1` or `1`.
//! [`ndarray`][ndarray]'s arrays are used to store spin lattice
//! and perform computations. The lattice can either be generated
//! by [`rand`][rand]s RNG or provided by the user.

pub mod lattice;
