use bitarr::BitArrNa;
use bitvec::prelude as bv;
use bitvec::slice::AsBits;
use dist_mat::DistMat;
use flate2::read::GzDecoder;
use float_cmp::approx_eq;
use itertools::{izip, Itertools};
use ndarray::{Array1, Array2};
use rayon::prelude::*;
use std::{error, fmt, fs, io, io::prelude::*};

mod bitarr;
pub mod cli;
mod dist_mat;

#[derive(fmt::Debug, Clone)]
pub enum GtWeights {
    None,
    All(Option<Array2<f32>>),
    Single(Option<Array2<f32>>),
}

#[derive(fmt::Debug, Clone)]
pub enum StrainsWith {
    G1,
    P1G1,
}

#[derive(fmt::Debug)]
pub struct Args {
    pub snps_fname: String,
    pub phen_fname: String,
    pub dists_fname: String,
    pub max_iter: usize,
    pub log_every: usize,
    pub delta: f64,
    pub gt_weights: GtWeights,
    pub rel_gt_weight: f32,
    pub phen_weights: Option<Array2<f32>>,
    pub rel_phen_weight: f32,
    pub p0g1_extra_weight: f32,
    pub dist_weight: f32,
    pub p1g1_filter: u32,
    pub avg_dist_strains: StrainsWith,
    pub threads: usize,
    pub out_fname: Option<String>,
    pub log_fname: Option<String>,
    pub avg_dists_fname: Option<String>,
    pub args_string: String,
}

#[derive(fmt::Debug)]
pub struct Calc {
    pub gt_weights: GtWeights,
    pub rel_gt_weight: f32,
    pub phen_weights: Option<Array2<f32>>,
    pub rel_phen_weight: f32,
    pub p0g1_extra_weight: f32,
    pub avg_dist_strains: StrainsWith,
    pub dist_weight: f32,
    pub p1g1_filter: u32,
    pub snps_arr: Vec<BitArrNa>,
    pub snps_var_ids: Vec<String>,
    pub phen: bv::BitVec,
    pub dists: DistMat<f32>,
    pub avg_dists: Option<Array1<f32>>,
    pub scores: Option<Array1<f64>>,
    pub active: Option<Vec<usize>>,
}

impl Calc {
    pub fn from_args(args: &Args) -> Calc {
        let dists = DistMat::from_csv_symmetric(&args.dists_fname).unwrap_or_else(|err| {
            eprintln!(
                "Error parsing distance file \"{}\": {}",
                args.dists_fname, err
            );
            std::process::exit(1);
        });
        let (phen, phen_sample_ids) = read_phen(&args.phen_fname).unwrap_or_else(|err| {
            eprintln!(
                "Error parsing phenotype file \"{}\": {}",
                args.phen_fname, err
            );
            std::process::exit(1);
        });
        let (snps_sample_ids, snps_var_ids, snps_arr) =
            read_snps(&args.snps_fname).unwrap_or_else(|err| {
                eprintln!("Error parsing SNPs file \"{}\": {}", args.snps_fname, err);
                std::process::exit(1);
            });
        // make sure that all 3 input files featured the same samples
        if !(phen_sample_ids == dists.labels && dists.labels == snps_sample_ids) {
            eprintln!(
                "Error: Sample labels in the distance file \"{}\", \
                genotype file \"{}\", and phenotype file \"{}\" are not equal.",
                &args.dists_fname, &args.snps_fname, &args.phen_fname
            );
            std::process::exit(1);
        }
        Calc {
            gt_weights: args.gt_weights.clone(),
            rel_gt_weight: args.rel_gt_weight,
            phen_weights: args.phen_weights.clone(),
            rel_phen_weight: args.rel_phen_weight,
            p0g1_extra_weight: args.p0g1_extra_weight,
            avg_dist_strains: args.avg_dist_strains.clone(),
            dist_weight: args.dist_weight,
            p1g1_filter: args.p1g1_filter,
            snps_arr,
            snps_var_ids,
            phen,
            dists,
            avg_dists: None,
            scores: None,
            active: None,
        }
    }

    /// Get the average pairwise p1g1 distance for every SNP in `snps_arr`.
    /// Then, normalize by mean.
    pub fn get_relative_avg_dists(&mut self) {
        let avg_dists: Vec<f32> = self
            .snps_arr
            .par_iter()
            .map(|snps| self.get_avg_dist(&snps))
            .collect();

        // normalize by all pairwise inter-sample distances:
        // let overall_mean = self.dists.mean() as f32;
        // or instead, the original way of normalizing:
        let overall_mean = avg_dists.iter().sum::<f32>() / avg_dists.len() as f32;

        self.avg_dists = Some(avg_dists.iter().map(|x| x / overall_mean).collect());
    }

    /// Get average pairwise distance among the relevant strains (i.e. with the respective
    /// genotype and phenotype or only the genotype, depending on `self.avg_dist_strains`).
    fn get_avg_dist(&self, snps: &BitArrNa) -> f32 {
        let width = get_bitvec_width(&snps.bits) as usize;
        let mut indices: Vec<usize> = Vec::new();
        match self.avg_dist_strains {
            // as closures can currently not be forced to be inlined, best performance
            // in some situations can only be achieved by either duplicating loops
            // (as below) or defining extra helper functions and inlining those.
            StrainsWith::P1G1 => {
                for (i, (g1, p1)) in izip!(snps.bits.as_slice(), self.phen.as_slice()).enumerate() {
                    let p1g1 = p1 & g1;
                    if p1g1 != 0 {
                        for (j, bit) in p1g1.bits::<bv::Local>().iter().enumerate() {
                            if *bit {
                                indices.push(i * width + j);
                            }
                        }
                    }
                }
            }
            StrainsWith::G1 => {
                for (i, g1) in snps.bits.as_slice().iter().enumerate() {
                    if *g1 != 0 {
                        for (j, bit) in g1.bits::<bv::Local>().iter().enumerate() {
                            if *bit {
                                indices.push(i * width + j);
                            }
                        }
                    }
                }
            }
        }
        self.dists.avg_pairwise_dist(&indices) as f32
    }

    pub fn write_avg_dists(&self, fname: &str) -> Result<(), Box<dyn error::Error>> {
        let dists = self
            .avg_dists
            .as_ref()
            .expect("Error: No average distances calculated yet.");
        // open file
        let mut file = fs::File::create(fname)?;
        writeln!(file, "varID,avg_pw_dist")?;
        for (dist, var_id) in dists.iter().zip(&self.snps_var_ids) {
            writeln!(file, "{},{}", var_id, dist)?;
        }
        Ok(())
    }

    pub fn prepare_weights(&mut self) {
        if self.phen_weights.is_none() {
            let mut w0 = self.phen.len() as f32 / 2f32 / self.phen.count_zeros() as f32;
            let mut w1 = self.phen.len() as f32 / 2f32 / self.phen.count_ones() as f32;
            if !approx_eq!(f32, self.rel_phen_weight, 1f32, ulps = 2) {
                w0 = w0 * self.rel_phen_weight + (1f32 - self.rel_phen_weight);
                w1 = w1 * self.rel_phen_weight + (1f32 - self.rel_phen_weight);
            }
            self.phen_weights = Some(ndarray::arr2(&[[w0], [w1]]));
            // self.phen_weights = Some((w0, w1));
        }
        match self.gt_weights {
            GtWeights::Single(None) => {
                let mut weights = Array2::zeros((2, self.snps_arr.len()));
                self.snps_arr.iter().enumerate().for_each(|(i, snps)| {
                    let zeros = snps.count_zeros();
                    let ones = snps.count_ones();
                    let sum = zeros + ones;
                    let w0 = sum as f32 / 2f32 / zeros as f32;
                    let w1 = sum as f32 / 2f32 / ones as f32;
                    weights[[0, i]] = w0;
                    weights[[1, i]] = w1;
                });

                if !approx_eq!(f32, self.rel_gt_weight, 1f32, ulps = 2) {
                    weights = weights * self.rel_gt_weight + (1f32 - self.rel_gt_weight);
                }

                self.gt_weights = GtWeights::Single(Some(weights));
            }
            GtWeights::All(None) => {
                let (mut zeros, mut ones) = (0, 0);
                self.snps_arr.iter().for_each(|snps| {
                    zeros += snps.count_zeros();
                    ones += snps.count_ones();
                });

                let sum = zeros + ones;
                let w0 = sum as f32 / 2f32 / zeros as f32;
                let w1 = sum as f32 / 2f32 / ones as f32;
                let mut weights = ndarray::arr2(&[[w0], [w1]]);

                if !approx_eq!(f32, self.rel_gt_weight, 1f32, ulps = 2) {
                    weights = weights * self.rel_gt_weight + (1f32 - self.rel_gt_weight);
                }
                self.gt_weights = GtWeights::All(Some(weights));
            }
            _ => (),
        }
    }

    pub fn get_scores(&mut self) {
        let phen_weights = self
            .phen_weights
            .as_ref()
            .expect("Error: Cannot calculate scores with unset phen_weights.");
        let gt_weights = match &self.gt_weights {
            GtWeights::Single(Some(weights)) => weights,
            GtWeights::All(Some(weights)) => weights,
            _ => panic!(
                "Error: Cannot calculate scores with unset gt_weights: {:?}.",
                self.gt_weights
            ),
        };
        let avg_pw_dists = self
            .avg_dists
            .as_ref()
            .expect("Error: Cannot calculate scores with unset avg_dists.");

        let (p1g1s, _, _, p0g1s) = self.get_all_pg_counts();
        let mut counts =
            &p1g1s * &phen_weights.slice(ndarray::s![1, ..]) * gt_weights.slice(ndarray::s![1, ..])
                - &p0g1s
                    * &phen_weights.slice(ndarray::s![0, ..])
                    * gt_weights.slice(ndarray::s![1, ..])
                    * self.p0g1_extra_weight;

        // filter counts to remove SNPs with fewer p1g1 than `p1g1_filter`
        ndarray::par_azip!((count in &mut counts, &p1g1 in &p1g1s){
            if (p1g1 < self.p1g1_filter as f32) || (*count < 0.)  {
                *count = 0.;
            }
        });

        let mut scores: Array1<f64> = Array1::zeros(counts.dim());
        ndarray::par_azip!((score in &mut scores, &count in &counts, &dist in avg_pw_dists)
            *score = (count + (dist - 1.) * count * self.dist_weight).into());

        self.scores = Some(scores);
    }

    pub fn get_all_pg_counts(&self) -> (Array1<f32>, Array1<f32>, Array1<f32>, Array1<f32>) {
        let mut p1g1s = Vec::new();
        let mut p0g0s = Vec::new();
        let mut p1g0s = Vec::new();
        let mut p0g1s = Vec::new();
        for snps in &self.snps_arr {
            let (p1g1, p0g0, p1g0, p0g1) = self.pg_counts(&snps);
            p1g1s.push(p1g1);
            p0g0s.push(p0g0);
            p1g0s.push(p1g0);
            p0g1s.push(p0g1);
        }
        (
            Array1::from(p1g1s),
            Array1::from(p0g0s),
            Array1::from(p1g0s),
            Array1::from(p0g1s),
        )
    }

    fn pg_counts(&self, snps: &BitArrNa) -> (f32, f32, f32, f32) {
        let mut p1g1 = 0f32;
        let mut p0g0 = 0f32;
        let mut p1g0 = 0f32;
        let mut p0g1 = 0f32;

        for (g1, nn, p1) in izip!(
            snps.bits.as_slice(),
            snps.not_nas.as_slice(),
            self.phen.as_slice()
        ) {
            let g0 = !g1 & nn;
            let p0 = !p1;
            p1g1 += (p1 & g1).count_ones() as f32;
            p0g0 += (p0 & g0).count_ones() as f32;
            p1g0 += (p1 & g0).count_ones() as f32;
            p0g1 += (p0 & g1).count_ones() as f32;
        }
        (p1g1, p0g0, p1g0, p0g1)
    }

    pub fn hhs_update_scores(
        &mut self,
        delta: f64,
        max_iter: usize,
        log_every: usize,
        logfile_fname: Option<&String>,
    ) {
        // initiliaze log file if there will be logging
        let mut logfile = None;
        if let Some(fname) = logfile_fname {
            if log_every > 0 {
                logfile = Some(fs::File::create(&fname).unwrap_or_else(|err| {
                    eprintln!("Error opening log file \"{}\": {}", fname, err);
                    std::process::exit(1);
                }));
            }
        }

        let orig_scores = self
            .scores
            .as_ref()
            .expect("Error: Cannot run HHS iteration with unset scores.")
            .clone();

        let mut scores: Vec<f64> = Vec::with_capacity(orig_scores.len());
        let mut active: Vec<usize> = Vec::with_capacity(orig_scores.len());

        for (i, &score) in orig_scores.iter().enumerate() {
            if score > 0.0 {
                scores.push(score);
                active.push(i);
            }
        }
        let sum: f64 = scores.iter().sum();
        // we will set negative scores to 0.0 later and subsequently drop scores that
        // are <= 0.0. in order to not trip over float-comparisons, we will use this
        // threshold instead of zero:
        let zero_threshold = sum / scores.len() as f64 * 1e-6;

        let mut new_scores: Vec<f64> = vec![0.; scores.len()];
        scores.shrink_to_fit();
        active.shrink_to_fit();

        let mut n_iter: usize = 0;
        loop {
            // write progress to STDOUT
            if n_iter % 10000 == 0 {
                println!(
                    "{} iterations done -- SNPs remaining: {}",
                    n_iter,
                    scores.len()
                );
            }
            // potential logging
            if let Some(file) = logfile.as_mut() {
                if n_iter % log_every == 0 {
                    write!(file, "{} iterations:", n_iter).expect("Error writing output");
                    log_scores(file, &scores, &self.snps_var_ids, &active);
                }
            }
            new_scores
                .par_iter_mut()
                .zip(&active)
                .zip(&scores)
                .for_each(|((new_score, &idx_i), current_score)| {
                    // the p1g1g1 with itself is the p1g1
                    let p1g1 = self.get_p1g1g1(&self.snps_arr[idx_i], &self.snps_arr[idx_i]);
                    // scale the score by the prevalence of the genotype--phenotype interaction
                    let mut tmp_score = current_score * p1g1;
                    // the score decreases according to all the overlaps with the other variants
                    for (&idx_j, other_score) in active.iter().zip(scores.iter()) {
                        if idx_i == idx_j {
                            continue;
                        }
                        tmp_score -= other_score
                            * self.get_p1g1g1(&self.snps_arr[idx_i], &self.snps_arr[idx_j])
                            * delta;
                    }
                    // no negative scores allowed
                    if tmp_score < 0.0 {
                        tmp_score = 0.0;
                    }
                    *new_score = tmp_score / p1g1;
                });
            // remove the indices of SNPs whose score fell below 0 from `active`
            for i in (0..active.len()).rev() {
                // to not get tripped by float-comparisons just use a small value instead of zero
                if new_scores[i] <= zero_threshold {
                    active.swap_remove(i);
                    new_scores.swap_remove(i);
                }
            }
            // rescale the new scores
            let norm_factor = sum / new_scores.iter().sum::<f64>();
            new_scores.iter_mut().for_each(|x| *x *= norm_factor);

            // every 1,000th iterations compare the old and new scores and stop if they
            // didn't change (i.e. convergence) or the maximum number of iterations (as
            // specified by the user) was reached.
            if n_iter % 1000 == 0 && scores.len() == new_scores.len() {
                // get sum of absolute differences
                let diffsum = scores
                    .iter()
                    .zip(&new_scores)
                    .fold(0.0, |sum, (s1, s2)| sum + (s1 - s2).abs());
                if diffsum < zero_threshold || n_iter >= max_iter {
                    break;
                }
            }
            // update new scores before next iteration
            scores.truncate(new_scores.len());
            scores.copy_from_slice(&new_scores);
            n_iter += 1;
        }

        if n_iter >= max_iter {
            println!(
                "\nMax. number of iterations ({}) exceeded: {} SNPs remain\n",
                max_iter,
                scores.len()
            );
        } else {
            println!(
                "\nConvergence after {}-{} iterations: {} SNPs remain\n",
                n_iter - 1000,
                n_iter,
                scores.len()
            );
        }

        // if there's logging, also log the final result
        if let Some(file) = logfile.as_mut() {
            write!(file, "{} iterations: ", n_iter).expect("Error writing output");
            log_scores(file, &scores, &self.snps_var_ids, &active);
        }

        // due to the use of swap_remove above, the SNPs are no longer sorted --> sort again
        let mut zipped = active.iter().zip(scores).collect::<Vec<_>>();
        zipped.sort_by_key(|(i, _)| *i);
        let (active, scores): (Vec<usize>, Vec<f64>) = zipped.into_iter().unzip();
        self.scores = Some(Array1::from(scores));
        self.active = Some(active);
    }

    fn get_p1g1g1(&self, snps_a: &BitArrNa, snps_b: &BitArrNa) -> f64 {
        let mut p1g1g1 = 0f64;
        for (p1, g1_a, g1_b) in izip!(
            self.phen.as_slice(),
            snps_a.bits.as_slice(),
            snps_b.bits.as_slice(),
        ) {
            p1g1g1 += (p1 & g1_a & g1_b).count_ones() as f64;
        }
        p1g1g1
    }

    /// Write scores to a buffer in csv format
    pub fn write_scores_csv<Buffer: io::Write>(&self, buffer: &mut Buffer) {
        let scores = self
            .scores
            .as_ref()
            .expect("Error: No scores calculated yet.");
        let indices = self
            .active
            .as_ref()
            .expect("Error: No scores calculated yet.");
        let dists = self
            .avg_dists
            .as_ref()
            .expect("Error: No scores calculated yet.");
        // check if there are scores to write and write them
        writeln!(buffer, "VarID,score,p1g1,p0g1,dist").expect("Error writing output");
        if (scores.len() == indices.len()) && !scores.is_empty() {
            for (i, score) in indices.iter().zip(scores) {
                let (p1g1, _, _, p0g1) = self.pg_counts(&self.snps_arr[*i]);
                let dist = dists[*i];
                writeln!(
                    buffer,
                    "{},{:.2},{},{},{:.3}",
                    self.snps_var_ids[*i], score, p1g1, p0g1, dist
                )
                .unwrap_or_else(|_| panic!("Error writing score with index {}", i));
            }
        }
    }
}

/// Write scores to a buffer in a single line in `index: score` format
pub fn log_scores<Buffer, T>(
    buffer: &mut Buffer,
    scores: &[T],
    var_ids: &[String],
    indices: &[usize],
) where
    Buffer: io::Write,
    T: fmt::Display,
{
    // check if there are scores to write and write them
    if (scores.len() == indices.len()) & !scores.is_empty() {
        write!(
            buffer,
            "{}",
            indices
                .iter()
                .zip(scores)
                .map(|(i, score)| format!("{}:{}", var_ids[*i], score))
                .collect::<Vec<String>>()
                .join(",")
        )
        .expect("Error logging score");
    }
    writeln!(buffer).expect("Error logging score");
}

/// read input file with snps in rows and samples in columns like
/// 010X11001
/// 01100X010
/// where anything other than '0' or '1' specifies unknown values
pub fn read_snps(
    infname: &str,
) -> Result<(Vec<String>, Vec<String>, Vec<BitArrNa>), Box<dyn error::Error>> {
    let mut snps_arr: Vec<BitArrNa> = Vec::new();
    // box file object to deal with different types from
    // `fs::File::open()` and `GzDecoder::new()`
    let mut file: Box<dyn io::Read> = Box::new(fs::File::open(infname)?);
    if infname.ends_with(".gz") {
        file = Box::new(GzDecoder::new(file))
    }
    let mut lines = io::BufReader::new(file).lines();
    // the first line holds the comma-separated sample names
    let first_line = lines.next().unwrap_or_else(|| {
        eprintln!("Error: Looks like the SNP file \"{}\" is empty", infname);
        std::process::exit(1)
    })?;
    let samples: Vec<String> = first_line
        .trim()
        .split(',')
        .map(|smp| smp.to_string())
        .collect();
    let num_samples = samples.len();
    let mut var_ids: Vec<String> = Vec::new();
    // now iterate over the remaining lines. they have the format "VarID:0010101010X01..."
    for (i, line) in lines.enumerate() {
        if let Ok(line) = line {
            // split the line into the variant label and the string with the actual genotypes
            let (var_id, snps_str) =
                line.trim()
                    .splitn(2, ':')
                    .collect_tuple()
                    .unwrap_or_else(|| {
                        eprintln!(
                            "Error: Looks like the format in the {}th line of the SNP file \
                    \"{}\" did not match what was expected (e.g. \"VarID:0010101001...\"). \
                    The first 50 characters were \"{}\".",
                            i + 2,
                            infname,
                            &line[..50]
                        );
                        std::process::exit(1)
                    });
            var_ids.push(var_id.to_string());
            // check whether there are as many genotypes as there were samples in the first line
            if snps_str.len() != num_samples {
                eprintln!(
                    "Error parsing genotype file \"{}\": The number of genotypes for the \
                variant \"{}\" was different from the number of samples provided in the \
                first line ({} vs. {}).",
                    infname,
                    var_id,
                    snps_str.len(),
                    num_samples
                );
                std::process::exit(1);
            }
            let snp = BitArrNa::from_string(&snps_str);
            snps_arr.push(snp);
        } else {
            eprintln!(
                "Error reading SNP input file \"{}\" at line {}.",
                infname,
                i + 2
            );
            std::process::exit(1);
        }
    }
    Ok((samples, var_ids, snps_arr))
}

/// reads csv with phenotype data. agnostic to whether there is a header.
fn read_phen(infname: &str) -> Result<(bv::BitVec, Vec<String>), Box<dyn error::Error>> {
    // box file object to deal with different types from
    // `fs::File::open()` and `GzDecoder::new()`
    let mut file: Box<dyn io::Read> = Box::new(fs::File::open(infname)?);
    if infname.ends_with("gz") {
        file = Box::new(GzDecoder::new(file))
    }
    // assume that there is a header but check if the second field in the header
    // is '0' or '1' in case there is none.
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_reader(file);
    // initiliase vector of strings for the labels and a string for the
    // phenotypes, which is used to generate the bitvec later
    let mut labels: Vec<String> = Vec::new();
    let mut phen_str = String::new();
    let header = reader.headers()?;
    if &header[1] == "0" || &header[1] == "1" {
        // there was no header --> add the data to the labels and string
        labels.push(header[0].to_string());
        phen_str.push_str(&header[1].to_string());
    }
    for rec in reader.records() {
        let row = rec?;
        labels.push(row[0].to_string());
        let phen = &row[1];
        if !(phen == "0" || phen == "1") {
            // the phenotype file should only hold '0' or '1' in the second column
            eprintln!(
                "Error parsing phenotype file \"{}\": The phenotype \
                for \"{}\" was \"{}\", but should be '0' or '1'.",
                infname, &row[0], phen
            );
            std::process::exit(1);
        }
        phen_str.push_str(phen);
    }
    // now generate bitvec and return
    let mut phen_bv = bv::bitvec![0; phen_str.len()];
    for (i, c) in phen_str.chars().enumerate() {
        // we would have thrown an error if there was something other than '0' or
        // '1' in `phen_str`. therefore, we only need to check for '1' and set
        // the bits accordingly.
        if c == '1' {
            phen_bv.set(i, true);
        }
    }
    Ok((phen_bv, labels))
}

/// Get number of bits in a `BitVec`'s storage unit.
fn get_bitvec_width<O, T>(_: &bv::BitVec<O, T>) -> u8
where
    O: bitvec::order::BitOrder,
    T: bitvec::store::BitStore,
{
    T::BITS
}

pub fn max(arr: &[f32]) -> f32 {
    arr.iter().fold(0f32, |a, &b| a.max(b))
}

pub fn max_a(arr: &Array1<f32>) -> f32 {
    arr.iter().fold(0f32, |a, &b| a.max(b))
}

pub fn max2(arr: &Array2<f32>) -> f32 {
    arr.iter().fold(0f32, |a, &b| a.max(b))
}
