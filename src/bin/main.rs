use hhs::*;

fn main() {
    let threads = 4;
    // create thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .expect("Error initializing threadpool");

    let args = Args::default();
    let mut calc = Calc::from_args(&args);
    calc.prepare_weights();
    calc.get_p1g1_relative_avg_dists();
    calc.get_scores();
    calc.hhs_update_scores(1e-4, 1e4 as usize);
    // calc.hhs_update_scores_old(1e-4, 1e4 as usize);
}
