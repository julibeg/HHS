use crate::{Args, GtWeights, StrainsWith};

pub fn parse_cmd_line() -> Args {
    let matches = clap::App::new("Hungry, Hungry SNPos")
        .about(
            "Run HHS on a binary genotype matrix, pairwise inter-sample distances, \
            and binary phenotypes (e.g. antibiotic resistance). \
            The algorithm finds variants associated with the phenotype while removing \
            false-positives caused by overlapping genotypes (e.g. co-occurrent resistance).\n\n\

            The input file holding the variants should look like:\n\
            ---------------------------------------\n\
            SampleID1,SampleID2,SampleID3,...\n\
            VarID1:10X0101X...\n\
            VarID2:X00X1001...\n\
            ...\n\
            ---------------------------------------\n\
            with the sample names in the first line and the genotypes of the single variants \
            below. Any character other than '0' or '1' can represent missing values.\n\n\
            
            The phenotypes should be provided in a CSV file with '0' and '1' denoting the \
            absence or presence of the phenotype, respectively. Missing values are not \
            allowed in the phenotype file and such samples should be removed from all \
            files before running the program.\n\n\

            To account for population structure, a symmetric pairwise distance matrix of the \
            samples should be provided in a CSV file with the sample names in the header (i.e. \
            the first row) and in the index column (and '0's in the main diagonal).
            ",
        )
        .version(clap::crate_version!())
        .arg(
            clap::Arg::with_name("GTs")
                .help("genotypes input file")
                .takes_value(true)
                .short("g")
                .long("gt")
                .required(true)
                .value_name("FILE")
                .display_order(1),
        )
        .arg(
            clap::Arg::with_name("phen")
                .help("phenotype input file")
                .takes_value(true)
                .short("p")
                .long("phen")
                .required(true)
                .value_name("FILE")
                .display_order(2),
        )
        .arg(
            clap::Arg::with_name("dist")
                .help("pairwise distances input file")
                .takes_value(true)
                .short("d")
                .long("dist")
                .required(true)
                .value_name("FILE")
                .display_order(3),
        )
        .arg(
            clap::Arg::with_name("max_iter")
                .help("maximum number of iterations")
                .takes_value(true)
                .short("i")
                .long("iter")
                .value_name("NUM")
                .default_value("3e5")
                .display_order(6),
        )
        .arg(
            clap::Arg::with_name("delta")
                .help("effect size multiplier (higher --> faster, but coarser results)")
                .takes_value(true)
                .short("D")
                .long("delta")
                .value_name("NUM")
                .default_value("1e-3")
                .display_order(7),
        )
        .arg(
            clap::Arg::with_name("output")
                .help("output file; if missing, result is printed to STDOUT.")
                .takes_value(true)
                .short("o")
                .long("output")
                .value_name("FILE")
                .display_order(8),
        )
        .arg(
            clap::Arg::with_name("threads")
                .help("number of threads; '0' will use all available CPUs")
                .takes_value(true)
                .short("t")
                .long("threads")
                .default_value("1")
                .value_name("NUM")
                .display_order(9),
        )
        .arg(
            clap::Arg::with_name("gt_weights")
                .help("how to calculate genotype weights")
                .takes_value(true)
                .long("gt_weights")
                .value_name("STR")
                .possible_values(&["all", "single"])
                .default_value("single")
                .display_order(10),
        )
        .arg(
            clap::Arg::with_name("rel_gt_weight")
                .help("relative genotype weight")
                .takes_value(true)
                .long("rel_gt_weights")
                .value_name("NUM")
                .default_value("1")
                .display_order(11),
        )
        // .arg(
        //     clap::Arg::with_name("phen_weights")
        //         .help("phenotype weights")
        //         .takes_value(true)
        //         .long("phen_weights")
        //         .value_name("STR")
        //         .default_value("1/1"),
        // )
        .arg(
            clap::Arg::with_name("rel_phen_weight")
                .help("relative phenotpye weight")
                .takes_value(true)
                .long("rel_phen_weights")
                .value_name("NUM")
                .default_value("1")
                .display_order(12),
        )
        .arg(
            clap::Arg::with_name("p0g1_extra_weight")
                .help("p0g1 extra weight")
                .takes_value(true)
                .long("p0g1_extra_weight")
                .value_name("NUM")
                .default_value("1")
                .display_order(13),
        )
        .arg(
            clap::Arg::with_name("dist_weight")
                .help("distance weight")
                .takes_value(true)
                .long("dist_weight")
                .value_name("NUM")
                .default_value("1")
                .display_order(14),
        )
        .arg(
            clap::Arg::with_name("p1g1_filter")
                .help("p1g1 filter")
                .takes_value(true)
                .long("p1g1_filter")
                .value_name("NUM")
                .default_value("2")
                .display_order(15),
        )
        .arg(
            clap::Arg::with_name("p1g1_avg_dist")
                .help("use avg. pairwise distances among all p1g1 strains instead of all g1 strains")
                .long("p1g1_avg_dist")
                .display_order(16),
        )
        .arg(
            clap::Arg::with_name("log_every")
                .help(
                    "write intermediary results every NUM iterations to the output file. \
                If '0', no logging will take place. Otherwise, a log file must be specified.",
                )
                .takes_value(true)
                .long("log_every")
                .value_name("NUM")
                .default_value("0")
                .display_order(17),
        )
        .arg(
            clap::Arg::with_name("logfile")
                .help("log file; required if `log_every` is other than '0'")
                .takes_value(true)
                .long("logfile")
                .value_name("FILE")
                .display_order(18),
        )
        .arg(
            clap::Arg::with_name("avg_dists_file")
                .help("write the avg. pairwise distances to this file before iterative elimination")
                .takes_value(true)
                .long("avg_dists_file")
                .value_name("FILE")
                .display_order(19),
        )
        .get_matches();

    // it's safe to call `unwrap` on arguments `required` by clap.
    let snps_fname = matches.value_of("GTs").unwrap().to_string();
    let phen_fname = matches.value_of("phen").unwrap().to_string();
    let dists_fname = matches.value_of("dist").unwrap().to_string();
    let max_iter =
        clap::value_t!(matches.value_of("max_iter"), f32).unwrap_or_else(|e| e.exit()) as usize;
    let log_every =
        clap::value_t!(matches.value_of("log_every"), f32).unwrap_or_else(|e| e.exit()) as usize;
    let delta = clap::value_t!(matches.value_of("delta"), f32).unwrap_or_else(|e| e.exit());
    let gt_weights = match matches.value_of("gt_weights").unwrap() {
        "all" => GtWeights::All(None),
        "single" => GtWeights::Single(None),
        _ => {
            // redundant error handling required for avoiding type mismatch. The validity of the
            // arguments has already been checked by clapped.
            eprintln!("Argument for 'gt_weights' must be 'all' or 'single'");
            std::process::exit(1);
        }
    };
    let rel_gt_weight =
        clap::value_t!(matches.value_of("rel_gt_weight"), f32).unwrap_or_else(|e| e.exit());
    let rel_phen_weight =
        clap::value_t!(matches.value_of("rel_phen_weight"), f32).unwrap_or_else(|e| e.exit());
    let p0g1_extra_weight =
        clap::value_t!(matches.value_of("p0g1_extra_weight"), f32).unwrap_or_else(|e| e.exit());
    let dist_weight =
        clap::value_t!(matches.value_of("dist_weight"), f32).unwrap_or_else(|e| e.exit());
    let p1g1_filter =
        clap::value_t!(matches.value_of("p1g1_filter"), u32).unwrap_or_else(|e| e.exit());
    let threads = clap::value_t!(matches.value_of("threads"), usize).unwrap_or_else(|e| e.exit());
    let avg_dist_strains = if matches.is_present("p1g1_avg_dist") {
        StrainsWith::P1G1
    } else {
        StrainsWith::G1
    };
    let out_fname = match matches.value_of("output") {
        Some(fname) => Some(fname.to_string()),
        None => None,
    };
    let log_fname = match matches.value_of("logfile") {
        Some(fname) => Some(fname.to_string()),
        None => {
            // when log_every != 0 we need a log file --> throw an error if there is none
            if log_every == 0 {
                None
            } else {
                eprintln!(
                    "Error: A log file is required (You specified that intermediary \
                    results should be logged every {} iterations.)",
                    log_every
                );
                std::process::exit(1);
            }
        }
    };
    let avg_dists_fname = match matches.value_of("avg_dists_file") {
        Some(fname) => Some(fname.to_string()),
        None => None,
    };
    // get command line args as string
    let args_string: String = std::env::args().collect::<Vec<String>>().join(" ");

    Args {
        snps_fname,
        phen_fname,
        dists_fname,
        max_iter,
        log_every,
        delta,
        gt_weights,
        rel_gt_weight,
        phen_weights: None,
        rel_phen_weight,
        p0g1_extra_weight,
        dist_weight,
        p1g1_filter,
        avg_dist_strains,
        threads,
        out_fname,
        log_fname,
        avg_dists_fname,
        args_string,
    }
}
