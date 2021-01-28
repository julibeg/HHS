use bitarr::BitArrNa;
use bitvec::prelude as bv;
use bitvec::slice::AsBits;
use dist_mat::DistMat;
use float_cmp::approx_eq;
use itertools::izip;
use ndarray::{Array1, Array2};
use rayon::prelude::*;
use std::fmt;
use std::fs;
use std::io;
use std::io::prelude::*;

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
    pub snps_na_char: char,
    pub snps_file_transposed: bool,
    pub phen_fname: String,
    pub dists_fname: String,
    pub iterations: usize,
    pub log_every: usize,
    pub delta: f32,
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
    pub phen: bv::BitVec,
    pub dists: DistMat<f32>,
    pub avg_dists: Option<Array1<f32>>,
    pub scores: Option<Array1<f32>>,
    pub active: Option<Vec<usize>>,
}

impl Calc {
    pub fn from_args(args: &Args) -> Calc {
        let dists = DistMat::from_csv_symmetric(&args.dists_fname).unwrap_or_else(|err| {
            eprintln!(
                "Error processing distance file '{}': {}",
                args.dists_fname, err
            );
            std::process::exit(1);
        });
        Calc {
            gt_weights: args.gt_weights.clone(),
            rel_gt_weight: args.rel_gt_weight,
            phen_weights: args.phen_weights.clone(),
            rel_phen_weight: args.rel_phen_weight,
            p0g1_extra_weight: args.p0g1_extra_weight,
            avg_dist_strains: args.avg_dist_strains.clone(),
            dist_weight: args.dist_weight,
            p1g1_filter: args.p1g1_filter,
            snps_arr: if args.snps_file_transposed {
                read_snps_transposed(&args.snps_fname, args.snps_na_char)
            } else {
                read_snps(&args.snps_fname, args.snps_na_char)
            },
            phen: read_phen(&args.phen_fname),
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

    pub fn write_avg_dists(&self, fname: &str) {
        let dists = self
            .avg_dists
            .as_ref()
            .expect("Error: No average distances calculated yet.");
        // check if there are scores to write and write them
        let mut file = fs::File::create(fname).unwrap_or_else(|err| {
            eprintln!(
                "Error opening file '{}' to write avg. pairwise distances to: {}",
                fname, err
            );
            std::process::exit(1);
        });
        for dist in dists.iter() {
            writeln!(file, "{}", dist).unwrap_or_else(|err| {
                panic!("Error writing avg. pairwise distance '{}': {}", dist, err)
            });
        }
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

        let mut scores: Array1<f32> = Array1::zeros(counts.dim());
        ndarray::par_azip!((score in &mut scores, &count in &counts, &dist in avg_pw_dists)
            *score = count + (dist - 1.) * count * self.dist_weight);

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
        delta: f32,
        n_iter: usize,
        log_every: usize,
        logfile_fname: Option<&String>,
    ) {
        // initiliaze log file if there will be logging
        let mut logfile = None;
        if let Some(fname) = logfile_fname {
            if log_every > 0 {
                logfile = Some(fs::File::create(&fname).unwrap_or_else(|err| {
                    eprintln!("Error opening log file '{}': {}", fname, err);
                    std::process::exit(1);
                }));
            }
        }

        let orig_scores = self
            .scores
            .as_ref()
            .expect("Error: Cannot run HHS iteration with unset scores.")
            .clone();

        let mut scores: Vec<f32> = Vec::with_capacity(orig_scores.len());
        let mut active: Vec<usize> = Vec::with_capacity(orig_scores.len());

        for (i, &score) in orig_scores.iter().enumerate() {
            if score > 0. {
                scores.push(score);
                active.push(i);
            }
        }
        let sum: f32 = scores.iter().sum();

        let mut new_scores: Vec<f32> = vec![0.; scores.len()];
        scores.shrink_to_fit();
        active.shrink_to_fit();

        for n in 0..n_iter {
            // write progress to STDOUT
            if n % (n_iter / 10) == 0 {
                println!("{} iterations done -- SNPs remaining: {}", n, scores.len());
            }

            // potential logging
            if let Some(file) = logfile.as_mut() {
                if n % log_every == 0 {
                    write!(file, "{} iterations: ", n).expect("Error writing output");
                    write_scores_single_line(file, &scores, &active);
                }
            }
            new_scores
                .par_iter_mut()
                .zip(&active)
                .zip(&scores)
                .for_each(|((new_score, &idx_i), current_score)| {
                    let p1g1 = self.get_p1g1g1(&self.snps_arr[idx_i], &self.snps_arr[idx_i]);

                    let mut tmp_score = current_score * p1g1;
                    for (&idx_j, other_score) in active.iter().zip(scores.iter()) {
                        if idx_i == idx_j {
                            continue;
                        }
                        tmp_score -= other_score
                            * self.get_p1g1g1(&self.snps_arr[idx_i], &self.snps_arr[idx_j])
                            * delta;
                    }
                    if tmp_score < 0. {
                        tmp_score = 0.;
                    }
                    *new_score = tmp_score / p1g1;
                });

            // remove the indices of SNPs whose score fell below 0 from `active`
            for i in (0..active.len()).rev() {
                if new_scores[i] <= 1e-5 {
                    active.swap_remove(i);
                    new_scores.swap_remove(i);
                }
            }

            // update and rescale scores before next iteration
            scores.truncate(new_scores.len());
            scores.copy_from_slice(&new_scores);
            let norm_factor = sum / scores.iter().sum::<f32>();
            scores.iter_mut().for_each(|x| *x *= norm_factor);
        }

        println!(
            "\nAll {} iterations done: {} SNPs remain\n",
            n_iter,
            scores.len()
        );

        // if there's logging, also log the final result
        if let Some(file) = logfile.as_mut() {
            write!(file, "{} iterations: ", n_iter).expect("Error writing output");
            write_scores_single_line(file, &scores, &active);
        }

        // due to the use of swap_remove above, the SNPs are no longer sorted --> sort again
        let mut zipped = active.iter().zip(scores).collect::<Vec<_>>();
        zipped.sort_by_key(|(i, _)| *i);
        let (active, scores): (Vec<usize>, Vec<f32>) = zipped.into_iter().unzip();
        self.scores = Some(Array1::from(scores));
        self.active = Some(active);
    }

    fn get_p1g1g1(&self, snps_a: &BitArrNa, snps_b: &BitArrNa) -> f32 {
        let mut p1g1g1 = 0f32;
        for (p1, g1_a, g1_b) in izip!(
            self.phen.as_slice(),
            snps_a.bits.as_slice(),
            snps_b.bits.as_slice(),
        ) {
            p1g1g1 += (p1 & g1_a & g1_b).count_ones() as f32;
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
        writeln!(buffer, "SNP_idx, score, p1g1, p0g1, dist").expect("Error writing output");
        if (scores.len() == indices.len()) & !scores.is_empty() {
            for (i, score) in indices.iter().zip(scores) {
                let (p1g1, _, _, p0g1) = self.pg_counts(&self.snps_arr[*i]);
                let dist = dists[*i];
                writeln!(
                    buffer,
                    "{}, {:.2}, {}, {}, {:.3}",
                    i, score, p1g1, p0g1, dist
                )
                .unwrap_or_else(|_| panic!("Error writing score with index {}", i));
            }
        }
    }
}

/// Write scores to a buffer in a single line in `index: score` format
pub fn write_scores_single_line<Buffer, T1, T2>(buffer: &mut Buffer, scores: &[T1], indices: &[T2])
where
    Buffer: io::Write,
    T1: fmt::Display,
    T2: fmt::Display,
{
    // check if there are scores to write and write them
    if (scores.len() == indices.len()) & !scores.is_empty() {
        let mut zipped = indices.iter().zip(scores);
        let (i, score) = zipped.next().unwrap();
        write!(buffer, "{}: {}", i, score)
            .unwrap_or_else(|_| panic!("Error writing score with index {}", i));
        for (i, score) in zipped {
            write!(buffer, ", {}: {}", i, score)
                .unwrap_or_else(|_| panic!("Error writing score with index {}", i));
        }
    }
    writeln!(buffer).expect("Error writing output");
}

/// read input file with snps in rows and samples in columns like
/// 010X11001
/// 01100X010
/// where 'X' specifies unknown values
pub fn read_snps(infname: &str, na_char: char) -> Vec<BitArrNa> {
    let mut snps_arr: Vec<BitArrNa> = Vec::new();

    let infile = fs::File::open(&infname).unwrap_or_else(|err| {
        eprintln!("Error opening SNPs input file '{}': {}", infname, err);
        std::process::exit(1);
    });
    let lines = io::BufReader::new(infile).lines();

    for (i, line) in lines.enumerate() {
        if let Ok(line) = line {
            let snp = BitArrNa::from_string(&line, na_char).unwrap_or_else(|err| {
                eprintln!(
                    "Error generating SNP-bitarr at input file '{}' line {}: {}",
                    infname,
                    i + 1,
                    err
                );
                std::process::exit(1);
            });
            snps_arr.push(snp);
        } else {
            eprintln!("Error reading input file '{}' at line {}", infname, i + 1);
            std::process::exit(1);
        }
    }
    snps_arr
}

/// read input file with snps in columns and samples in rows like
/// 00
/// 11
/// 01
/// X0
/// 10
/// 1X
/// 00
/// 01
/// 10
/// where 'X' specifies unknown values
pub fn read_snps_transposed(infname: &str, na_char: char) -> Vec<BitArrNa> {
    // read file first to determine dimensions
    let infile = fs::File::open(&infname).unwrap_or_else(|err| {
        eprintln!("Error opening SNPs input file '{}': {}", infname, err);
        std::process::exit(1);
    });
    let mut lines = io::BufReader::new(infile).lines();
    // get first line and count number of SNPs
    let n_snps = match lines.next() {
        Some(Ok(line)) => line.len(),
        _ => {
            eprintln!("Error parsing SNPs input file. Is it empty?");
            std::process::exit(1);
        }
    };
    // get number of lines --> number of samples
    let n_samples = lines.count() + 1; // add 1 for the first line

    // initialize empty snp_arr with `BitArrNa`s of the correct length
    let mut snps_arr: Vec<BitArrNa> = (0..n_snps).map(|_| BitArrNa::new(n_samples)).collect();

    // now the buffer's iterator has been consumed and we have to open the file again. we can call
    // `unwrap` here, because it has already been opened without error.
    let infile = fs::File::open(&infname).unwrap();
    let lines = io::BufReader::new(infile).lines();

    for (i, line) in lines.enumerate() {
        if let Ok(sample) = line {
            for (j, c) in sample.chars().enumerate() {
                if c == '0' {
                    continue;
                } else if c == '1' {
                    snps_arr[j].bits.set(i, true);
                } else if c == na_char {
                    snps_arr[j].not_nas.set(i, false);
                } else {
                    eprintln!(
                        "Error parsing SNPs file '{}': Char at position {} was '{}'; expected '0', '1' or '{}'.",
                        infname,
                        j + 1,
                        c,
                        na_char
                    );
                    std::process::exit(1);
                }
            }
        } else {
            eprintln!("Error reading input file '{}' at line {}", infname, i + 1);
            std::process::exit(1);
        }
    }
    snps_arr
}

pub fn read_phen(infname: &str) -> bv::BitVec {
    let phen_str = fs::read_to_string(infname).unwrap_or_else(|err| {
        eprintln!("Error opening input file '{}': {}", infname, err);
        std::process::exit(1);
    });
    let phen_str = phen_str.trim();
    let mut phen = bv::bitvec![0; phen_str.len()];

    for (i, c) in phen_str.chars().enumerate() {
        if c == '0' {
            continue;
        } else if c == '1' {
            phen.set(i, true);
        } else {
            eprintln!(
                "Error parsing phenotype file '{}': Char at position {} was '{}'; \
                     expected '0' or '1'",
                infname,
                i + 1,
                c
            );
            std::process::exit(1);
        }
    }
    phen
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
