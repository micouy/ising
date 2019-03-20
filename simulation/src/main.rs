#![allow(non_snake_case)]

use ::chrono::prelude::*;
use ::ndarray::prelude::*;
use ::pbr::ProgressBar;
use ::rand::prelude::*;
use ::rayon::prelude::*;
use ::serde_json::{json, to_string_pretty};

use std::{
    env::args,
    fs,
    path::Path,
    sync::mpsc::{channel, Sender},
    thread,
};

use ising_lib::prelude::*;

const SIZE: usize = 50;
const T_MIN: f64 = 0.1;
const T_MAX: f64 = 5.0;
const T_STEP: f64 = 0.1;
const FLIPS_TO_SKIP: usize = 60_000;
const MEASUREMENTS_PER_T: usize = 1000;
const ATTEMPTS_PER_FLIP: usize = 20;

struct Params {
    T_range: (f64, f64),
    flips_to_skip: usize,
    measurements_per_T: usize,
    flips_per_measurement: usize,
    attempts_per_flip: usize,
    lattice_size: usize,
    J: f64,
    h: f64,
}

impl Params {
    fn new(J: f64, h: f64) -> Self {
        Self {
            T_range: (T_MIN, T_MAX),
            flips_to_skip: FLIPS_TO_SKIP,
            measurements_per_T: MEASUREMENTS_PER_T,
            flips_per_measurement: SIZE * SIZE,
            attempts_per_flip: ATTEMPTS_PER_FLIP,
            lattice_size: SIZE,
            J,
            h,
        }
    }
}

struct Record {
    T: f64,
    dE: f64,
    I: f64,
    X: f64,
}

fn compose_results(records: &[Record], params: Params) -> String {
    let records = records
        .iter()
        .map(|r| {
            json!({
                "T": r.T,
                "I": r.I,
                "dE": r.dE,
                "X": r.X,
            })
        })
        .collect::<Vec<_>>();

    to_string_pretty(&json!({
        "records": records,
        "params": {
            "J": params.J,
            "h": params.h,
        },
    }))
    .unwrap()
}

fn compose_file_name() -> String {
    let now = Local::now().format("%d.%m.%Y-%H.%M").to_string();
    let id = thread_rng().gen_range(100_i32, 999_i32);

    format!("results-{}-{}.txt", now, id)
}

fn cmp_by_T(a: &Record, b: &Record) -> std::cmp::Ordering {
    a.T.partial_cmp(&b.T).unwrap_or(std::cmp::Ordering::Less)
}

fn run(params: Params, pb_tx: Sender<()>) -> (String, String) {
    let mut rng = SmallRng::from_entropy();
    let mut lattice = Lattice::new((params.lattice_size, params.lattice_size));
    let Ts: Vec<f64> = TRange::from_step(params.T_range.0, params.T_range.1, T_STEP).collect();
    let h = Array::from_elem((params.lattice_size, params.lattice_size), params.h);

    // "cool" the lattice to its natural state
    (0..params.flips_to_skip).for_each(|_| {
        let _ = (0..params.attempts_per_flip)
            .map(|_| {
                let ix = lattice.gen_random_index();
                let E_diff = lattice.measure_E_diff(ix, params.J);
                let probability = calc_flip_probability(E_diff, params.T_range.0);

                if probability > rng.gen() {
                    lattice.flip_spin(ix);

                    true
                } else {
                    false
                }
            })
            .take_while(|already_flipped| !already_flipped)
            .count();
    });

    let mut records: Vec<Record> = Ts
        .into_iter()
        .map(|T| {
            let (Es, Is) = (0..params.measurements_per_T)
                .map(|_| {
                    (0..params.flips_per_measurement).for_each(|_| {
                        let _ = (0..params.attempts_per_flip)
                            .map(|_| {
                                let ix = lattice.gen_random_index();
                                let E_diff = lattice.measure_E_diff_with_h(ix, &h, params.J);
                                let probability = calc_flip_probability(E_diff, T);

                                if probability > rng.gen() {
                                    lattice.flip_spin(ix);

                                    true // the flip has already occured
                                } else {
                                    false // the flip has not occured yet
                                }
                            })
                            .take_while(|already_flipped| !already_flipped)
                            .count();
                    });

                    let _ = pb_tx.send(());

                    (lattice.measure_E(params.J), lattice.measure_I())
                })
                .unzip::<_, _, Vec<_>, Vec<_>>();

            let dE = calc_dE(&Es, T);
            let I = calc_I(&Is);
            let X = calc_X(&Es);

            Record { T, dE, I, X }
        })
        .collect();

    let file_name = compose_file_name();
    records.sort_by(cmp_by_T);
    let results = compose_results(&records, params);

    (results, file_name)
}

fn main() {
    let dir_name = args()
        .nth(1)
        .expect("Specify the directory you want to save the results to.");

    // make sure it's a valid directory
    assert!(Path::new(&dir_name).is_dir());

    let Js = vec![0.2, 0.6, 1.0, 1.4, 1.8];

    let hs = vec![0.4, 0.8, 1.2, 1.6, 2.0];

    let Js_and_hs = Js
        .into_iter()
        .map(|J| hs.clone().into_iter().map(move |h| (J, h)))
        .flatten()
        .collect::<Vec<_>>();

    let bar_count =
        ((T_MAX - T_MIN) / T_STEP).floor() as u64 * MEASUREMENTS_PER_T as u64 * Js_and_hs.len() as u64;
    let (pb_tx, pb_rx) = channel();

    let handle = thread::spawn(move || {
        let mut pb = ProgressBar::new(bar_count);
        pb.set_width(Some(100));
        pb.show_message = true;
        pb.message("Running...");

        for _ in 0..bar_count {
            let _ = pb_rx.recv();
            pb.inc();
        }

        pb.finish_print("Finished!");
    });

    let results = Js_and_hs
        .into_iter()
        .zip((0..).map(|_| pb_tx.clone()))
        .collect::<Vec<_>>()
        .into_par_iter()
        .map(|((J, h), pb_tx)| {
            let params = Params::new(J, h);

            run(params, pb_tx)
        })
        .collect::<Vec<(String, String)>>();

    let _ = handle.join();

    for (result, file_name) in results {
        let path = format!("{}/{}", dir_name, file_name);
        fs::write(path, result.as_bytes()).unwrap();
    }
}
