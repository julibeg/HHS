use hhs::*;

fn main() {
    let args = cli::parse_cmd_line();

    // create thread pool
    let threads = args.threads;
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .expect("Error initializing threadpool");

    // println!("{:?}", args);

    let mut calc = Calc::from_args(&args);

    // for snp in calc.snps_arr {
    //     println!("{}\n-------", snp);
    // }
    // println!("{:?}", calc.snps_arr);

    calc.prepare_weights();
    calc.get_p1g1_relative_avg_dists();
    calc.get_scores();
    calc.hhs_update_scores(1e-4, 1e6 as usize);
}
