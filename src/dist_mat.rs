use std::cmp;
use std::fs;
use std::io;
use std::ops;

/// Triangular distance matrix
#[derive(Debug)]
pub struct DistMat<T> {
    pub mat: Vec<Vec<T>>,
}

impl<T> DistMat<T> {
    /// initialize empty DistMat with given capacity
    pub fn new(n: usize) -> DistMat<T> {
        let mut mat: Vec<Vec<T>> = Vec::with_capacity(n);
        for i in 0..n {
            mat.push(Vec::with_capacity(n - i));
        }
        DistMat { mat }
    }
}

impl<T> DistMat<T>
where
    T: num::Num + std::string::ToString + Copy,
{
    /// Write DistMat to buffer in square/symmetric form
    pub fn write_csv_symmetric<Buffer: io::Write>(&self, buffer: &mut Buffer) {
        let n = self.mat.len();
        for i in 0..=n {
            let mut line: Vec<T> = Vec::with_capacity(n);
            for j in 0..=n {
                let dist: T;
                if i == j {
                    dist = T::zero();
                } else {
                    dist = self[[i, j]];
                }
                line.push(dist);
            }
            let line: Vec<String> = line.into_iter().map(|i| i.to_string()).collect();
            writeln!(buffer, "{}", &line.join(","))
                .unwrap_or_else(|_| panic!("Error writing result at i: {}", i));
        }
    }
}

impl<T> DistMat<T>
where
    T: num::Num,
{
    /// generate DistMat from .csv file
    pub fn from_csv_symmetric(infname: &str) -> Result<DistMat<T>, io::Error> {
        let infile = fs::File::open(infname).unwrap_or_else(|err| {
            eprintln!("Error opening input file {}: {}", infname, err);
            std::process::exit(1);
        });
        let infile = io::BufReader::new(infile);
        DistMat::read_symmetric(infile, ",")
    }

    /// parse square/symmetric input into DistMat
    pub fn read_symmetric<Buffer: io::BufRead>(
        buffer: Buffer,
        delim: &str,
    ) -> Result<DistMat<T>, io::Error> {
        let mut size: Option<usize> = None; // is set after reading first line
        let mut mat: Vec<Vec<T>> = Vec::new();
        for (i, line) in buffer.lines().enumerate() {
            let line = line?;
            let mut vec = Vec::new();
            if let Some(size) = size {
                // check if size has already been set
                vec.reserve(size - i); // and allocate space if so
            }
            for elem in line.split(delim).skip(i + 1) {
                let elem = match T::from_str_radix(elem.trim(), 10) {
                    Ok(elem) => elem,
                    Err(_) => {
                        // from_str_radix returns a fmt error
                        return Err(io::Error::new(
                            // --> a new io::Error needs to be created
                            io::ErrorKind::InvalidInput,
                            format!(
                                "Error converting input \"{}\" into number at line {}",
                                elem, i
                            ),
                        ));
                    }
                };
                vec.push(elem);
            }
            if i == 0 {
                // we know the size after reading the first
                size = Some(vec.len()) // line --> can be used to pre-allocate
            }; // memory for the following lines

            mat.push(vec);
        }
        mat.truncate(mat.len() - 1); // the last vec was empty due to skip
        Ok(DistMat { mat }) // remove it from the matrix
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
