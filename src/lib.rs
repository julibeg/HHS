use crate::bitarr::BitArrNa;
use crate::dist_mat::DistMat;
use bitvec::prelude as bv;
use itertools::izip;
use std::fs;
use std::io;
use std::io::BufRead;

pub mod bitarr;
pub mod dist_mat;

enum gt_weights {
    snp,
    all,
    None,
    Some((f32, f32))
}

pub struct args {
    gt_weights: gt_weights,
    rel_gt_weight: Option<f32>,
    phen_weights: Option<(f32, f32)>,
    rel_phen_weight: Option<f32>,
    P0G1_extra_weight: f32,
    dist_weight: f32,
    p1g1_filter: u32,
}

/// read input file with snps in rows and samples in columns like
/// 010X11001
/// 01100X010
/// and 'X' specifying unknown values
pub fn read_snps(infname: &str, na_char: char) -> Vec<BitArrNa> {
    let mut all_snps: Vec<BitArrNa> = Vec::new();

    let infile = fs::File::open(&infname).unwrap_or_else(|err| {
        eprintln!("Error opening input file {}: {}", infname, err);
        std::process::exit(1);
    });
    let infile = io::BufReader::new(infile);

    for (i, line) in infile.lines().enumerate() {
        if let Ok(line) = line {
            let snp = BitArrNa::from_string(&line, na_char).unwrap_or_else(|err| {
                eprintln!(
                    "Error generating SNP-bitarr at input file {} line {}: {}",
                    infname,
                    i + 1,
                    err
                );
                std::process::exit(1);
            });
            all_snps.push(snp);
        } else {
            eprintln!("Error reading input file at line {}", i + 1);
            std::process::exit(1);
        }
    }
    all_snps
}

pub fn pg_counts(snps: &BitArrNa, phen: &bv::BitVec) -> (u64, u64, u64, u64) {
    let mut p1g1 = 0u64;
    let mut p0g0 = 0u64;
    let mut p1g0 = 0u64;
    let mut p0g1 = 0u64;

    for (g1, nn, p1) in izip!(
        snps.bits.as_slice(),
        snps.not_nas.as_slice(),
        phen.as_slice()
    ) {
        let g0 = !g1 & nn;
        let p0 = !p1;
        p1g1 += (p1 & g1).count_ones() as u64;
        p0g0 += (p0 & g0).count_ones() as u64;
        p1g0 += (p1 & g0).count_ones() as u64;
        p0g1 += (p0 & g1).count_ones() as u64;
    }
    (p1g1, p0g0, p1g0, p0g1)
}

pub fn read_phen(infname: &str) -> bv::BitVec {
    let phen_str = fs::read_to_string(infname).unwrap_or_else(|err| {
        eprintln!("Error opening input file {}: {}", infname, err);
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
                "Error parsing phenotype file {}: Char at position {} was \'{}\'; \
                 expected \'0\' or \'1\'",
                infname,
                i + 1,
                c
            );
            std::process::exit(1);
        }
    }
    phen
}
