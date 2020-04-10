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
    T: num_traits::Num + std::string::ToString + Copy,
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
                    // sort i and j to figure out where to look up the value
                    let smaller = cmp::min(i, j);
                    let larger = cmp::max(i, j);
                    dist = self[smaller][larger - smaller - 1];
                }
                line.push(dist);
            }
            let line: Vec<String> = line.into_iter().map(|i| i.to_string()).collect();
            writeln!(buffer, "{}", &line.join(","))
                .expect(&format!("Error writing result at i: {}", i));
        }
    }
}

impl<T> DistMat<T>
where
    T: num_traits::Num,
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
    T: num_traits::Num + Copy + std::fmt::Debug + num_traits::cast::ToPrimitive,
{
    /// get the average pairwise distance among the samples specified by `indices`.
    /// assumes `indices` to be sorted. goes over every index in indices and considers all
    /// distances between the corresponding sample and the other samples represented by the
    /// remaining indices.
    /// it does this by extracting the `vec` at the index's position in `mat`. `vec` holds the
    /// distances of the sample at index to all downstream samples (i.e. with higher indices).
    /// the matching distances are extracted by converting the indeces of the samples (i.e. the
    /// actual indices) into indices in `vec` by subtracting the index of the original sample
    /// and 1.
    pub fn avg_pairwise_dist(&self, indices: &[usize]) -> f64 {
        let mut sum = T::zero();
        let mut count = T::zero();
        for (i, &idx_i) in (&indices[..indices.len() - 1]).iter().enumerate() {
            let vec = &self[idx_i];
            for j in &indices[i + 1..] {
                count = count + T::one();
                sum = sum + vec[j - idx_i - 1];
            }
        }
        let result = sum
            .to_f64()
            .expect("Error converting avg_pairwise_dist result to f64")
            / count
                .to_f64()
                .expect("Error converting avg_pairwise_dist result to f64");
        result
    }
}

impl<T> ops::Index<usize> for DistMat<T> {
    type Output = Vec<T>;

    fn index(&self, i: usize) -> &Self::Output {
        &self.mat[i]
    }
}

impl<T> ops::IndexMut<usize> for DistMat<T> {
    // IndexMut is a child of Index --> Output does not need to be set
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        &mut self.mat[i]
    }
}
