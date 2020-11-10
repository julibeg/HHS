use crate::{Args, GtWeights};

const ALPHANUM_EXCEPT_01: &[&str] = &[
    "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S",
    "T", "U", "V", "W", "X", "Y", "Z", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l",
    "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "2", "3", "4", "5", "6",
    "7", "8", "9",
];

pub fn parse_cmd_line() -> Args {
    let matches = clap::App::new("Hungry, Hungry SNPos")
        .about(
            "Run HHS on a list of SNPs, a pairwise distance matrix and binary phenotypes (e.g. \
            antibiotic resistance). \
            The algorithm finds SNPs associated with the positive phenotype while removing \
            false-positives caused by overlapping genotypes (e.g. co-occurrent resistance).\n\
            The input file holding the SNPs should look like: \n\n\
            10X0101X \n\
            X00X1001 \n\n\
            with one SNP per row (and one sample per column). 'X' denotes missing \
            values (any character other than '0' or '1' can be selected to represent NAs). \
            The file can also be transposed (one sample per row and one SNP per column) if the \
            -T flag is provided.\n\
            The input file with the phenotypes should have the same layout but feature only a \
            single line. Samples with missing phenotype values are not allowed (i.e. no 'X's in \
            the phenotype file) and should be removed from all input files prior analysis. \n\
            The symmetric pairwise distance matrix should be provided in a csv file with '0's in \
            the diagonal.
            ",
        )
        .version(clap::crate_version!())
        .arg(
            clap::Arg::with_name("SNPs")
                .help("SNPs input file")
                .takes_value(true)
                .short("s")
                .long("snps")
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
            clap::Arg::with_name("NA-char")
                .help("character [A-Za-z2-9] specifying missing values in SNPs file")
                .takes_value(true)
                .default_value("X")
                .possible_values(ALPHANUM_EXCEPT_01)
                .hide_possible_values(true)
                .short("n")
                .long("NA-value")
                .value_name("CHAR")
                .display_order(4),
        )
        .arg(
            clap::Arg::with_name("transposed")
                .help("use when SNPs input file is transposed (SNPs per column, samples per row)")
                .short("T")
                .long("transposed")
                .display_order(5),
        )
        .arg(
            clap::Arg::with_name("iterations")
                .help("number of iterations")
                .takes_value(true)
                .short("i")
                .long("iter")
                .value_name("NUM")
                .default_value("3e4")
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
                .help(
                    "output file; if missing, result is printed to STDOUT.",
                )
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
            clap::Arg::with_name("log_every")
                .help(
                    "Write intermediary results every NUM iterations to the output file. \
                If '0', no logging will take place. Otherwise, a log file must be specified.",
                )
                .takes_value(true)
                .long("log_every")
                .value_name("NUM")
                .default_value("0")
                .display_order(16),
        )
        .arg(
            clap::Arg::with_name("logfile")
                .help(
                    "log file; required if `log_every` is other than '0'",
                )
                .takes_value(true)
                .long("logfile")
                .value_name("FILE")
                .display_order(17),
        )
        .get_matches();

    // it's safe to call `unwrap` on arguments `required` by clap.
    let snps_fname = matches.value_of("SNPs").unwrap().to_string();
    let snps_na_char = matches.value_of("NA-char").unwrap().chars().next().unwrap();
    let snps_file_transposed = matches.is_present("transposed");
    let phen_fname = matches.value_of("phen").unwrap().to_string();
    let dists_fname = matches.value_of("dist").unwrap().to_string();
    let iterations =
        clap::value_t!(matches.value_of("iterations"), f32).unwrap_or_else(|e| e.exit()) as usize;
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
    // get command line args as string
    let args_string: String = std::env::args().collect::<Vec<String>>().join(" ");

    Args {
        snps_fname,
        snps_na_char,
        snps_file_transposed,
        phen_fname,
        dists_fname,
        iterations,
        log_every,
        delta,
        gt_weights,
        rel_gt_weight,
        phen_weights: None,
        rel_phen_weight,
        p0g1_extra_weight,
        dist_weight,
        p1g1_filter,
        threads,
        out_fname,
        log_fname,
        args_string,
    }
}
