use hhs::*;
use std::fs;

fn main() {
    let args = cli::parse_cmd_line();

    // create thread pool
    let threads = args.threads;
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .expect("Error initializing threadpool");

    // setup HHS calculator instance
    println!("Reading ipnut data...\n");
    let mut calc = Calc::from_args(&args);

    // prepare weights, distances and scores
    println!("Preparing weights, calculating average pairwise distances and filtering initial scores...\n");
    calc.prepare_weights();
    calc.get_p1g1_relative_avg_dists();
    calc.get_scores();

    // run HHS and convert result to string
    println!("Running HHS...");
    let active = calc.hhs_update_scores(args.delta, args.iterations as usize);
    let result_string = active
        .into_iter()
        .map(|i| i.to_string())
        .collect::<Vec<String>>()
        .join(", ");

    // print result to stdout or file
    match args.out_fname {
        Some(fname) => {
            println!("Writing result to {}\n", fname);
            fs::write(fname, result_string).unwrap_or_else(|err| {
                eprintln!("Error creating output file: {}", err);
                std::process::exit(1);
            });
        }
        None => println!("Indices of remaining SNPs:\n{}\n", result_string),
    };

    println!("Done");
}
