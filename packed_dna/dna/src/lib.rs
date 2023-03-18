//! A general-purpose genomics crate for dealing with DNA.
//!
//! Provide a representation of nucleotides sequence that is more memory efficient.
//! [PackedDna](crate::packed::PackedDna) are used to represent a dna sequence
//!
//! #Usage
//!
//! ```toml
//! [dependencies]
//! dna = { path = "../dna" }
//! ```
//!
//! Next:
//!
//! ```rust
//! use dna::{packed::PackedDna, Nuc};
//! use std::{str::FromStr};
//!
//! main() {
//!    // input a nucleotide sequence
//!     let dna = "ATCGGACCT";
//!
//!     // store the input sequence in PackedDna Struct
//!     let packed_dna = PackedDna::from_str(&dna);
//!     // return if the input is unsupported
//!     if packed_dna.is_err() {
//!         println!("unsupported input");
//!         return
//!     }
//!     let dna_vec= packed_dna.unwrap();
//!
//!     // get a particular nucleotide with index i
//!     let i = 3;
//!     let pi_nuc = dna_vec.get(i-1);
//!     if pi_nuc.is_err() {
//!         println!("fail to get the {}-th nucleotide",i);
//!         return
//!     }
//!     println!("The {}-th nucleotide is {:?}",i,pi_nuc.unwrap())
//! }
//! ```
//! 

#![warn(missing_docs)]
use std::{convert::TryFrom, fmt::Display, str::FromStr};

/// This module provides a PackedDna struct,
/// a case insensitive `FromStr` implementation,
/// a `FromIterator` implementation,
/// and a `fn get(&self, idx: usize) -> Nuc` getter for a particular nucleotide
pub mod packed {
    use super::*;
    use std::iter::FromIterator;

    /// 1. A representation that is more memory efficient that simply storing a vector of `Nuc`
    /// the nucleotide sequence is stored in a vector of u8 elements
    /// the length of the sequence is stored in a usize element
    #[derive(Debug, Clone, PartialEq, Eq, Hash)]
    pub struct PackedDna {
        nucs_vec: Vec<u8>,
        length: usize,
    }

    /// An error that can occur when parsing a nucleotide.
    #[derive(Debug, thiserror::Error)]
    #[error("failed to parse nucleotide from {0}")]
    pub struct ParsePackedNucsError<T: Display>(T);

    /// 2. A FromStr implementation
    impl FromStr for PackedDna {
        type Err = ParsePackedNucsError<char>;

        /// This function convert a string slice to a PackedDna struct
        /// # Input: a slice of string
        /// # Output: Result
        /// # Error
        /// when given string contains characters other than "A","T","C","G",
        /// the functions return ParsePackedNucsError
        fn from_str(s: &str) -> Result<Self, Self::Err> {
            // convert all characters to uppercase
            let upper = s.to_ascii_uppercase();
            // create a new empty u8 vector
            let mut packed_nucs: Vec<u8> = Vec::new();
            // create variable 'len' to store the length of the sequence
            let mut len: usize = 0usize;
            // create variable 'nucs_unit' to store the value for every 4 nucleotides into an u8 element
            let mut nucs_unit: u8 = 0u8;

            // range over all characters in the string slice
            for nuc in upper.chars() {
                // when a u8 element is filled with the value of 4 nucleotides, push the element into vector
                if len % 4 == 0 && len != 0 {
                    packed_nucs.push(nucs_unit);
                    nucs_unit = 0u8;
                }

                // each nucleotide is represented using 2 bits
                match nuc {
                    // 0b00 represents A
                    'A' => {
                        nucs_unit <<= 2u8;
                    }
                    // 0b01 represents C
                    'C' => {
                        nucs_unit = (nucs_unit << 2u8) + 1u8;
                    }
                    // 0b10 represents G
                    'G' => {
                        nucs_unit = (nucs_unit << 2u8) + 2u8;
                    }
                    // 0b11 represents T
                    'T' => {
                        nucs_unit = (nucs_unit << 2u8) + 3u8;
                    }
                    // other characters are not acceptable
                    _ => return Err(ParsePackedNucsError(nuc)),
                }

                len += 1;
            }

            // push the last element
            if len != 0 {
                packed_nucs.push(nucs_unit);
            }

            Ok(PackedDna {
                nucs_vec: packed_nucs,
                length: len,
            })
        }
    }

    // 3. A `FromIterator` implementation to construct it from an iterator over `Nuc`s
    impl FromIterator<Nuc> for PackedDna {
        /// This function convert a Vector of Nuc to a PackedDna struct
        /// # Input: a Vector of Nuc
        /// # Output: Result
        fn from_iter<I: IntoIterator<Item = Nuc>>(iter: I) -> Self {
            // create a new empty u8 vector
            let mut packed_nucs: Vec<u8> = Vec::new();
            // create variable 'len' to store the length of the sequence
            let mut len: usize = 0;
            // create variable 'nucs_unit' to store the value for every 4 nucleotides into an u8 element
            let mut nucs_unit: u8 = 0;

            // range over all Nucs in the vector
            for nuc in iter {
                if len % 4 == 0 && len != 0 {
                    packed_nucs.push(nucs_unit);
                    nucs_unit = 0u8;
                }

                match nuc {
                    Nuc::A => {
                        nucs_unit <<= 2u8;
                    }
                    Nuc::C => {
                        nucs_unit = (nucs_unit << 2u8) + 1u8;
                    }
                    Nuc::G => {
                        nucs_unit = (nucs_unit << 2u8) + 2u8;
                    }
                    Nuc::T => {
                        nucs_unit = (nucs_unit << 2u8) + 3u8;
                    }
                }

                len += 1;
            }

            if len != 0 {
                packed_nucs.push(nucs_unit);
            }

            PackedDna {
                nucs_vec: packed_nucs,
                length: len,
            }
        }
    }

    /// An error that can occur when getting for a particular nucleotide.
    #[derive(Debug, thiserror::Error)]
    #[error("failed to get nucleotide from {0}")]
    pub struct GetPackedNucsError(usize);

    /// 4. A `fn get(&self, idx: usize) -> Nuc` getter for a particular nucleotide
    impl PackedDna {
        /// This function get a particular nucleotide with index idx in the PackedDna struct
        /// # Input: Self(PackedDna struct), index:usize
        /// # Output: Result
        /// # Error
        /// when the given index in out of range,
        /// or the retrieved value cannot be converted into a Nuc type,
        /// the function return GetPackedNucsError
        pub fn get(&self, idx: usize) -> Result<Nuc, GetPackedNucsError> {
            if idx >= self.length {
                return Err(GetPackedNucsError(idx));
            }

            let max_full_unit = self.len() / 4;
            let element_idx = idx / 4;
            let bit_idx = idx % 4;
            let mut nuc = self.nucs_vec[element_idx];

            if idx >= max_full_unit * 4 {
                // idx in a unit that is not full
                nuc >>= 2 * (self.len() - idx - 1);
                nuc &= 3;
            } else {
                // idx in a unit that is full
                nuc >>= 2 * (3 - bit_idx);
                nuc &= 3;
            }

            match nuc {
                0u8 => Ok(Nuc::A),
                1u8 => Ok(Nuc::C),
                2u8 => Ok(Nuc::G),
                3u8 => Ok(Nuc::T),
                _ => Err(GetPackedNucsError(idx)),
            }
        }

        /// This function get the length of the sequence stored in PackedDna struct
        /// # Input: Self(PackedDna struct)
        /// # Output: length(usize)
        pub fn len(&self) -> usize {
            self.length
        }

        /// This function check if a struct is empty
        /// # Input: Self(PackedDna struct)
        /// # Output: bool
        pub fn is_empty(&self) -> bool {
            self.length != 0
        }

        /// This function get the length of u8 vector stored in PackedDna struct
        /// # Input: Self(PackedDna struct)
        /// # Output: length(usize)
        pub fn vec(&self) -> Vec<u8> {
            let mut vec_new: Vec<u8> = Vec::new();
            for i in &self.nucs_vec {
                vec_new.push(*i)
            }
            vec_new
        }
    }
}

/// A nucleotide
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Nuc {
    /// Adenine
    A,
    /// Cytosine
    C,
    /// Guanine
    G,
    /// Thymine
    T,
}

/// An error that can occur when parsing a nucleotide.
#[derive(Debug, thiserror::Error)]
#[error("failed to parse nucleotide from {0}")]
pub struct ParseNucError<T: Display>(T);

/// a Tryfrom implementation for Nuc
impl TryFrom<char> for Nuc {
    type Error = ParseNucError<char>;

    /// This function converts a character into a Nuc
    /// # Input: a character
    /// # Output: Result
    /// # Error
    /// when given character is not among ' A','T','C','G',
    /// the functions return ParseNucError
    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value.to_ascii_uppercase() {
            'A' => Ok(Self::A),
            'C' => Ok(Self::C),
            'G' => Ok(Self::G),
            'T' => Ok(Self::T),
            _ => Err(ParseNucError(value)),
        }
    }
}

/// a FromStr implementation for Nuc
impl FromStr for Nuc {
    type Err = ParseNucError<String>;

    /// This function convert a string slice to a PackedDna struct
    /// # Input: a slice of string
    /// # Output: Result
    /// # Error
    /// when given string is not among "A","T","C","G",
    /// the functions return ParseNucError
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let upper = s.to_ascii_uppercase();
        match upper.as_str() {
            "A" => Ok(Self::A),
            "C" => Ok(Self::C),
            "G" => Ok(Self::G),
            "T" => Ok(Self::T),
            _ => Err(ParseNucError(upper)),
        }
    }
}

/// This module test the functions mentioned above
#[cfg(test)]
mod tests {
    use super::*;
    use std::{iter::FromIterator, str::FromStr};
    #[test]
    fn tryfrom_char() {
        // test if each nucleotide can be converted to Nuc
        // test if the result is case insensitive
        let a = Nuc::A;
        assert!(a == Nuc::try_from('A').unwrap());
        assert!(a != Nuc::try_from('T').unwrap());
        assert!(a == Nuc::try_from('a').unwrap());

        let c = Nuc::C;
        assert!(c == Nuc::try_from('C').unwrap());
        assert!(c != Nuc::try_from('A').unwrap());
        assert!(c == Nuc::try_from('c').unwrap());

        let g = Nuc::G;
        assert!(g == Nuc::try_from('G').unwrap());
        assert!(g != Nuc::try_from('A').unwrap());
        assert!(g == Nuc::try_from('g').unwrap());

        let t = Nuc::T;
        assert!(t == Nuc::try_from('T').unwrap());
        assert!(t != Nuc::try_from('A').unwrap());
        assert!(t == Nuc::try_from('t').unwrap());

        // test illegal input
        assert!(Nuc::try_from('B').is_err());
        assert!(Nuc::try_from('b').is_err());
    }

    #[test]
    fn fromstr_nuc() {
        // test if each nucleotide string slice can be converted to Nuc
        // test if the result is case insensitive
        let a = Nuc::A;
        assert!(a == Nuc::from_str("A").unwrap());
        assert!(a != Nuc::from_str("T").unwrap());
        assert!(a == Nuc::from_str("a").unwrap());

        let c = Nuc::C;
        assert!(c == Nuc::from_str("C").unwrap());
        assert!(c != Nuc::from_str("A").unwrap());
        assert!(c == Nuc::from_str("c").unwrap());

        let g = Nuc::G;
        assert!(g == Nuc::from_str("G").unwrap());
        assert!(g != Nuc::from_str("A").unwrap());
        assert!(g == Nuc::from_str("g").unwrap());

        let t = Nuc::T;
        assert!(t == Nuc::from_str("T").unwrap());
        assert!(t != Nuc::from_str("A").unwrap());
        assert!(t == Nuc::from_str("t").unwrap());

        // test illegal input
        assert!(Nuc::from_str("B").is_err());
        assert!(Nuc::from_str("b").is_err());
        assert!(Nuc::from_str("ATC").is_err());
        assert!(Nuc::from_str("").is_err());
    }

    #[test]
    fn fromstr_packed_dna() {
        // test if the 4 nucleotides can be packed into an u8 element
        let case = packed::PackedDna::from_str("ATCG");
        assert!(case.is_ok());
        let unwrap_case = case.unwrap();
        assert!((&unwrap_case).vec() == vec![54u8]);
        assert!((&unwrap_case).len() == 4);

        // test if the storage is case-insensitive
        let case = packed::PackedDna::from_str("atcg");
        assert!(case.is_ok());
        let unwrap_case = case.unwrap();
        assert!((&unwrap_case).vec() == vec![54u8]);
        assert!((&unwrap_case).len() == 4);

        // test if more than 4 nucleotides can be packed into several u8 elements
        let case = packed::PackedDna::from_str("atcgGGACgact");
        assert!(case.is_ok());
        let unwrap_case = case.unwrap();
        assert!((&unwrap_case).vec() == vec![54u8, 161u8, 135u8]);
        assert!((&unwrap_case).len() == 12);

        // test illegal inputs
        let case = packed::PackedDna::from_str("ABCg");
        assert!(case.is_err());

        // test empty input
        let case = packed::PackedDna::from_str("");
        assert!(case.is_ok());
        let unwrap_case = case.unwrap();
        assert!((&unwrap_case).len() == 0);
    }

    #[test]
    fn fromiter_packed_dna() {
        let v = vec![Nuc::A, Nuc::T, Nuc::C, Nuc::G];
        let case = packed::PackedDna::from_iter(v);
        assert!(case.get(0).unwrap() == Nuc::A);
        assert!(case.get(1).unwrap() == Nuc::T);
        assert!(case.get(2).unwrap() == Nuc::C);
        assert!(case.get(3).unwrap() == Nuc::G);
    }

    #[test]
    fn get_packed_dna() {
        // test 4 nucleotides case
        let case = packed::PackedDna::from_str("ATGC").unwrap();
        assert!(case.get(0).unwrap() == Nuc::A);
        assert!(case.get(2).unwrap() == Nuc::G);
        // test index out of range
        assert!(case.get(4).is_err());

        // test more than 4 nucleotides case
        let case = packed::PackedDna::from_str("ATGCGGCTA").unwrap();
        assert!(case.get(4).unwrap() == Nuc::G);
        assert!(case.get(8).unwrap() == Nuc::A);
        // test index out of range
        assert!(case.get(9).is_err());

        // test empty input
        let case = packed::PackedDna::from_str("").unwrap();
        assert!(case.get(4).is_err());
    }
}
