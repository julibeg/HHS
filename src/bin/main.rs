use hhs::*;

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
    calc.hhs_update_scores(
        args.delta,
        args.iterations,
        args.log_every,
        args.log_fname.as_ref(),
    );

    // print result to stdout or file
    match args.out_fname {
        Some(fname) => {
            println!("Writing result to {}\n", fname);
            // first, write the command that was used to invoke the program
            std::fs::write(&fname, args.args_string + "\n").unwrap_or_else(|err| {
                eprintln!("Error creating output file: {}", err);
                std::process::exit(1);
            });
            // then, write scores, distances etc. of the final SNPs in CSV format
            let mut file = std::fs::OpenOptions::new().append(true).open(&fname).unwrap_or_else(|err| {
                eprintln!("Error opening output file: {}", err);
                std::process::exit(1);
            });
            calc.write_scores_csv(&mut file);
        }

        None => {
            println!("Writing result to STDOUT:\n");
            // first, write the command that was used to invoke the program
            println!("{}", args.args_string);
            // then, write scores, distances etc. of the final SNPs in CSV format
            calc.write_scores_csv(&mut std::io::stdout());
        }
    };

    println!("Done");
}
