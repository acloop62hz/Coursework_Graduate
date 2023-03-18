// a nucleotide counter
//
// If run with `nuccount --dna ACGTTT" it should print the following to stdout:
// ```
// Input: ACGTTT
//
// A: 1
// C: 1
// G: 1
// T: 3
// ```
//
// be sure to exit with informative error messages if the input is invalid

use dna::{packed::PackedDna, Nuc};
use std::str::FromStr;
use structopt::StructOpt;

/// Count the number of occurrences of each nucleotide in the provided DNA.
#[derive(Debug, StructOpt)]
struct Opts {
    /// The DNA sequence for which we should retrieve a nucleotide count.
    ///
    /// It is case insensitive but only nucleotides A, C, G and T are supported.
    #[structopt(short = "d", long, required = true)]
    dna: String,
}

fn main() {
    // get input
    let opts = Opts::from_args();
    let dna = opts.dna;
    println!("Input: {}", &dna);
    let input = PackedDna::from_str(&dna);
    // return if input contains illegal characters
    if input.is_err() {
        println!("unsupported input");
        return;
    }
    let packed_dna = input.unwrap();
    let length = packed_dna.len();

    // set counters for 4 types of nucleotides
    let mut a_count: usize = 0;
    let mut t_count: usize = 0;
    let mut c_count: usize = 0;
    let mut g_count: usize = 0;

    for i in 0..length {
        let nuc = packed_dna.get(i);
        if nuc.is_err() {
            println!("fail to get nucleotide with index {}", i);
            return;
        }
        match nuc.unwrap() {
            Nuc::A => a_count += 1,
            Nuc::C => c_count += 1,
            Nuc::G => g_count += 1,
            Nuc::T => t_count += 1,
        }
    }
    // print counts
    println!("A: {}", a_count);
    println!("C: {}", c_count);
    println!("G: {}", g_count);
    println!("T: {}", t_count);
}
