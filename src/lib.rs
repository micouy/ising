//! Everything you need to run
//! [Ising model](https://en.wikipedia.org/wiki/Ising_model) simulation.
//! Despite its simplicity, the simulation allows us to observe an interesting
//! physical phenomenon - phase transition. [`ndarray`][ndarray]'s arrays
//! are used to store the spin lattice and perform computations. Spin array
//! either be generated by [`rand`][rand]s RNG or provided by the user.
#![deny(missing_docs)]

pub mod calculations;
pub mod lattice;
pub mod prelude;
