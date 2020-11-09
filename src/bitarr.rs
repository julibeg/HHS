use bitvec::prelude as bv;
use itertools::izip;
use std::fmt;

/// Set trailing bits of a `BitVec` (i.e. bits in the last storage element not belonging to the
/// actual `BitVec`) to zero. Should be endian-agnostic.
fn zero_trailing_bits<O, T>(bitvec: &mut bv::BitVec<O, T>)
where
    O: bitvec::order::BitOrder,
    T: bitvec::store::BitStore,
{
    let width = T::BITS as usize;
    let numb_elements = bitvec.capacity() / width;
    bitvec.as_mut_slice()[numb_elements - 1] = T::FALSE;
    for i in bitvec.capacity() - width..bitvec.len() {
        bitvec.set(i, true);
    }
}

#[derive(fmt::Debug)]
pub struct BitArrNa {
    pub bits: bv::BitVec<bv::Local, usize>,
    pub not_nas: bv::BitVec<bv::Local, usize>,
}

impl BitArrNa {
    /// create new `BitArrNa` with all-zero bits and all-one not_nas
    pub fn new(size: usize) -> BitArrNa {
        let bits = bv::bitvec![0; size];
        let mut not_nas = bv::bitvec![1; size];

        // the trailing bits of the last storage element in `not_nas` are also '1's which might
        // lead to unexpected results --> set them to zero
        zero_trailing_bits(&mut not_nas);
        BitArrNa { bits, not_nas }
    }

    pub fn from_string(string: &str, na_char: char) -> Result<BitArrNa, String> {
        let mut bitarr = BitArrNa::new(string.len());
        for (i, c) in string.chars().enumerate() {
            if c == '0' {
                continue;
            } else if c == '1' {
                bitarr.bits.set(i, true);
            } else if c == na_char {
                bitarr.not_nas.set(i, false);
            } else {
                return Err(format!(
                    "Char at position {} was \'{}\'; expected \'0\', \'1\' or \'{}\'.",
                    i + 1,
                    c,
                    na_char
                ));
            }
        }

        Ok(bitarr)
    }

    pub fn dist<T>(&self, other: &BitArrNa) -> T
    where
        T: num::Num + num::cast::FromPrimitive,
    {
        let mut result: T = T::zero();

        for (bits, not_nas, other_bits, other_not_nas) in izip!(
            self.bits.as_slice(),
            self.not_nas.as_slice(),
            other.bits.as_slice(),
            other.not_nas.as_slice()
        ) {
            let res_bits = bits ^ other_bits;
            let res_not_nas = not_nas & other_not_nas;
            let incr = T::from_u32((res_bits & res_not_nas).count_ones())
                .expect("Error converting distance to requested type");
            result = result + incr;
        }

        result
    }

    pub fn len(&self) -> usize {
        self.bits.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn count_ones(&self) -> usize {
        self.bits.count_ones()
    }

    pub fn count_zeros(&self) -> usize {
        self.bits.count_zeros() - self.not_nas.count_zeros()
    }
}

impl fmt::Display for BitArrNa {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "bits:\t\t{}\nnot_nas:\t{}", self.bits, self.not_nas)
    }
}
