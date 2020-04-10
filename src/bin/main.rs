use bitvec::prelude as bv;
use hhs::dist_mat::DistMat;
use hhs::*;

fn main() {
    let dists_infname = "200samples_10000snps.dists";
    let snps_infname = "10000snps_200samples.gt";
    let phen_infname = "200samples.phen";

    let dists: DistMat<u32> = DistMat::from_csv_symmetric(dists_infname).unwrap_or_else(|err| {
        eprintln!("Error processing input file {}: {}", dists_infname, err);
        std::process::exit(1);
    });

    let snps = read_snps(snps_infname, '2');

    let phen = read_phen(phen_infname);

    let bits = bv::bitvec![0, 1, 0, 0, 1, 0, 0];
    let nns = bv::bitvec![0, 1, 1, 1, 1, 0, 1];
    let phen2 = bv::bitvec![1, 0, 0, 1, 1, 0, 0];
    let bla = bitarr::BitArrNa {
        bits: bits,
        not_nas: nns,
    };
    let res = pg_counts(&snps[4762], &phen);

    println!("{:?}", res);
}
