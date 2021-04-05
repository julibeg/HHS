use flate2::read::GzDecoder;
use std::{cmp, error, fs, io, ops, str};

/// Triangular distance matrix
#[derive(Debug)]
pub struct DistMat<T> {
    pub mat: Vec<Vec<T>>,
    pub labels: Vec<String>,
}

impl<T> DistMat<T> {
    /// initialize empty DistMat with given capacity
    pub fn new(n: usize) -> DistMat<T> {
        let mut mat: Vec<Vec<T>> = Vec::with_capacity(n);
        for i in 0..n {
            mat.push(Vec::with_capacity(n - i));
        }
        DistMat {
            mat,
            labels: Vec::with_capacity(n - 1),
        }
    }
}

impl<T> DistMat<T>
where
    T: num::Num + str::FromStr + std::fmt::Display,
    <T as std::str::FromStr>::Err: std::error::Error + 'static,
{
    pub fn from_csv_symmetric(infname: &str) -> Result<DistMat<T>, Box<dyn error::Error>> {
        // box file object to deal with different types from 
        // `fs::File::open()` and `GzDecoder::new()`
        let mut file: Box<dyn io::Read> = Box::new(fs::File::open(infname)?);
        if infname.ends_with("gz") {
            file = Box::new(GzDecoder::new(file))
        }
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .from_reader(file);
        // parse header and skip first element, which holds the label of the index column
        let labels: Vec<String> = reader
            .headers()?
            .iter()
            .skip(1)
            .map(|label| label.to_string())
            .collect();
        // the upper half of the symmetric matrix in the csv is stored as a
        // vector of vectors with decreasing lengths --> initialise `mat`
        let mut mat: Vec<Vec<T>> = Vec::new();
        // iterate over rows in csv
        for (i, record) in reader.records().enumerate() {
            // transform current row to iterator right away. as he reader only allows
            // rows of equal length anyway, we can unwrap when calling `.next()` later.
            let record = record?;
            let mut row = record.iter();
            // get the first element --> the row label
            let label = row.next().unwrap();
            // check if the label matches up with the corresponding field
            // in the header; throw error otherwise.
            if label != &labels[i] {
                eprintln!(
                    "Error: \"{}\" is expected to hold a symmetric matrix, but the\
                    header and row labels do not match up at position {} with \"{}\"\
                    in the header and \"{}\" as row label.",
                    infname, i, &label, &labels[i]
                );
                std::process::exit(1);
            }
            // skip elements before the main diagonal
            let mut row = row.skip(i);
            // check that the diagonal entry is zero; again we can unwrap because we know the
            // row must be long enough
            let diag_entry: T = row.next().unwrap().parse()?;
            if diag_entry != T::zero() {
                eprintln!(
                    "Error parsing \"{}\": The diagonal entry in the {}th row (with index \
                    \"{}\") should be zero but is \"{}\"",
                    infname,
                    i + 1,
                    label,
                    diag_entry
                );
                std::process::exit(1);
            }
            // Take the entries to the right of the diagonal and push them onto `mat`
            let dists = row
                .map(|dist| dist.parse::<T>())
                .collect::<Result<Vec<T>, _>>()?;
            // There are no entries to the right of the diagonal in the final row and
            // it is going to be empty --> don't append it to `mat`.
            if dists.len() > 0 {
                mat.push(dists);
            }
        }
        Ok(DistMat { mat, labels })
    }
}

impl<T> DistMat<T>
where
    T: num::Num + Copy + std::fmt::Debug + num::cast::ToPrimitive,
{
    /// Get the average pairwise distance among the samples specified by `indices`.
    /// assumes `indices` to be sorted. Goes over every index in indices and considers all
    /// distances between the corresponding sample and the other samples represented by the
    /// remaining indices.
    /// It does this by extracting the `vec` at the index's position in `mat`. `vec` holds the
    /// distances of the sample at every index to all downstream samples (i.e. with higher indices).
    /// The relevant distances are extracted by converting the indeces of the samples (i.e. the
    /// actual indices) into indices in `vec` by subtracting the index of the original sample
    /// and 1.
    pub fn avg_pairwise_dist(&self, indices: &[usize]) -> f64 {
        if indices.len() < 2 {
            // better solution would be to return an error for len == 0
            return 0.;
        }
        let mut sum = T::zero();
        let mut count = T::zero();
        for (i, &idx_i) in (&indices[..indices.len() - 1]).iter().enumerate() {
            let vec = &self[idx_i];
            for j in &indices[i + 1..] {
                count = count + T::one();
                sum = sum + vec[j - idx_i - 1];
            }
        }
        sum.to_f64().unwrap_or_else(|| {
            panic!(
                "Error converting avg_pairwise_dist result to f64; sum: {:?}",
                sum
            )
        }) / count.to_f64().unwrap_or_else(|| {
            panic!(
                "Error converting avg_pairwise_dist result to f64; count: {:?}",
                count
            )
        })
    }

    /// Get overall mean distance in matrix
    pub fn mean(&self) -> f64 {
        let mut sum = T::zero();
        let mut len = T::zero();
        for vec in &self.mat {
            for &elem in vec {
                len = len + T::one();
                sum = sum + elem;
            }
        }
        sum.to_f64().unwrap_or_else(|| {
            panic!(
                "Error converting avg_pairwise_dist result to f64; sum: {:?}",
                sum
            )
        }) / len.to_f64().unwrap_or_else(|| {
            panic!(
                "Error converting avg_pairwise_dist result to f64; len: {:?}",
                len
            )
        })
    }
}

impl<T> ops::Index<usize> for DistMat<T> {
    type Output = Vec<T>;

    #[inline]
    fn index(&self, i: usize) -> &Self::Output {
        &self.mat[i]
    }
}

impl<T> ops::IndexMut<usize> for DistMat<T> {
    // IndexMut is a child of Index --> Output does not need to be set

    #[inline]
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        &mut self.mat[i]
    }
}

impl<T> ops::Index<[usize; 2]> for DistMat<T> {
    type Output = T;

    /// Going to panic if a[0] == a[1].
    #[inline]
    fn index(&self, a: [usize; 2]) -> &Self::Output {
        let smaller = cmp::min(a[0], a[1]);
        let larger = cmp::min(a[0], a[1]);
        &self.mat[smaller][larger - smaller - 1]
    }
}

impl<T> ops::IndexMut<[usize; 2]> for DistMat<T> {
    // IndexMut is a child of Index --> Output does not need to be set

    /// Going to panic if a[0] == a[1].
    #[inline]
    fn index_mut(&mut self, a: [usize; 2]) -> &mut Self::Output {
        let smaller = cmp::min(a[0], a[1]);
        let larger = cmp::min(a[0], a[1]);
        &mut self.mat[smaller][larger - smaller - 1]
    }
}
