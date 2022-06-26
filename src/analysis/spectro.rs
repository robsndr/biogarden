use itertools::iproduct;
use lazy_static::lazy_static;
use std::collections::HashMap;

use crate::ds::sequence::Sequence;
use crate::ds::tile::Tile;
use crate::error::{BioError, Result};

// Define the spectral masses (in Units) of different amino-bases
lazy_static! {
    static ref monoisotopic_mass_table: HashMap<u8, f64> = HashMap::from([
        (b'F', 147.06841),
        (b'I', 113.08406),
        (b'V', 99.06841),
        (b'L', 113.08406),
        (b'S', 87.03203),
        (b'P', 97.05276),
        (b'M', 131.04049),
        (b'T', 101.04768),
        (b'A', 71.03711),
        (b'Y', 163.06333),
        (b'-', 0.0),
        (b'H', 137.05891),
        (b'N', 114.04293),
        (b'D', 115.02694),
        (b'Q', 128.05858),
        (b'K', 128.09496),
        (b'E', 129.04259),
        (b'C', 103.00919),
        (b'G', 57.02146),
        (b'R', 156.10111),
        (b'W', 186.07931)
    ]);
}

/// Calculates the monoisotopic mass of an amino acid
///
/// The monoisotopic mass is computed by using the principal isotope of each atom in a amino acid.
/// Its average mass is taken by taking the average mass of each atom in the molecule over all naturally appearing isotopes.
/// Monoisotopic mass is used more often than average mass in mass spectrometry.
/// The standard unit used in mass spectrometry for measuring mass is the atomic mass unit (dalton (Da))
///
/// # Arguments
/// * `protein` - protein sequence
///
/// # Example
/// ```
/// use biotech::analysis::spectro::weighted_mass;
/// use biotech::ds::sequence::Sequence;
///
/// let protein = Sequence::from("SKADYEK");
/// assert_eq!(weighted_mass(&protein).unwrap(), 821.3919199999999);
/// ```
pub fn weighted_mass(protein: &Sequence) -> Result<f64> {
    let mut mass: f64 = 0.0;
    for amino in protein {
        mass += monoisotopic_mass_table
            .get(amino)
            .ok_or(BioError::ItemNotFound)?;
    }
    Ok(mass)
}

/// Returns the frequency and value of the most frequently occurring shift between peaks of two mass spectrums
///
/// Mass spectrography is used to infer protein from their weight.
/// Comparing spectral representations can give useful insight into the similarity of corresponding particles.
/// A given protein will often occur nested in a more complex structure, which will shift the peaks of the mass spectrum.
/// In order to still be able to quantify the similarity of the spectra, the shift must be calculated.
///
/// # Arguments
/// * `s1`, `s2` - evaluated mass spectra
///
/// # Example
/// ```
/// use biotech::analysis::spectro::spectral_mass_shift;
///
/// let spec_1 = Vec::<f64>::from([
///     186.07931, 
///     287.12699, 
///     548.20532, 
///     580.18077, 
///     681.22845, 
///     706.27446, 
///     782.27613
/// ]);
/// 
/// let spec_2 = Vec::<f64>::from([
///     101.04768, 
///     158.06914, 
///     202.09536, 
///     318.09979, 
///     419.14747, 
///     463.17369
/// ]);
///
/// let (peak_count, shift) = spectral_mass_shift(&spec_1, &spec_2).unwrap();
///
/// assert_eq!(peak_count, 3);
/// assert_eq!(shift, 8504.0);
/// ```
pub fn spectral_mass_shift(s1: &[f64], s2: &[f64]) -> Result<(usize, f64)> {
    // Calculate minkovsky difference for s1 x s2
    // TODO: Introduce better approximation/rounding function
    let minkovsky_diff: Vec<i64> = iproduct!(s1.iter(), s2.iter())
        .map(|(v, x)| ((v - x) * 100.0).ceil() as i64)
        .collect();

    // Calculate mode of obtained set equvalent to mass shift
    let mut counts = HashMap::<i64, usize>::new();
    let max = minkovsky_diff
        .iter()
        .copied()
        .max_by_key(|&n| {
            let count = counts.entry(n).or_insert(0);
            *count += 1;
            *count
        })
        .unwrap_or(0);

    let max_freq = *counts.get(&max).ok_or(BioError::ItemNotFound)?;
    Ok((max_freq, max as f64))
}

/// Constructs the spectrum of weighted masses for all prefixes of a given sequence
///
/// # Arguments
/// * `seq` - sequence to generate prefix spectrum
///
/// # Example
/// ```
/// use biotech::analysis::spectro::prefix_spectrum;
/// use biotech::ds::sequence::Sequence;
///
/// let protein = Sequence::from("IASWMQS");
/// let spectrum = prefix_spectrum(&protein);
/// 
/// let result = Vec::<f64>::from([
///     113.08406, 
///     184.12117, 
///     271.1532, 
///     457.23251000000005, 
///     588.273, 
///     716.33158, 
///     803.36361
/// ]);
///
/// assert_eq!(spectrum.unwrap(), result);
/// ```
pub fn prefix_spectrum(seq: &Sequence) -> Result<Vec<f64>> {
    let mut weight_sum = 0_f64;
    let spectrum: Result<Vec<f64>> = seq
        .into_iter()
        .map(|x| {
            weight_sum += monoisotopic_mass_table
                .get(x)
                .ok_or(BioError::ItemNotFound)?;
            Ok(weight_sum)
        })
        .collect();
    spectrum
}

/// Constructs the spectrum of weighted masses for all suffixes of a given sequence
///
/// # Arguments
/// * `seq` - sequence to generate suffix spectrum
///
/// # Example
/// ```
/// use biotech::analysis::spectro::suffix_spectrum;
/// use biotech::ds::sequence::Sequence;
///
/// let protein = Sequence::from("IASWMQS");
/// let spectrum = suffix_spectrum(&protein);
/// 
/// let result = Vec::<f64>::from([
///     87.03203, 
///     215.09061000000003, 
///     346.13110000000006, 
///     532.21041, 
///     619.24244, 
///     690.27955, 
///     803.36361
/// ]);
///
/// assert_eq!(spectrum.unwrap(), result);
/// ```
pub fn suffix_spectrum(seq: &Sequence) -> Result<Vec<f64>> {
    let mut weight_sum = 0_f64;
    let spectrum: Result<Vec<f64>> = seq
        .into_iter()
        .rev()
        .map(|x| {
            weight_sum += monoisotopic_mass_table
                .get(x)
                .ok_or(BioError::ItemNotFound)?;
            Ok(weight_sum)
        })
        .collect();
    spectrum
}

/// Constructs the complete mass spectrum for a given sequence (prefix + suffix spectrum)
///
/// # Arguments
/// * `seq` - sequence to generate complete spectrum
///
/// # Example
/// ```
/// use biotech::analysis::spectro::complete_spectrum;
/// use biotech::ds::sequence::Sequence;
///
/// let protein = Sequence::from("IASW");
/// let spectrum = complete_spectrum(&protein);
/// let result = Vec::<f64>::from([
///     113.08406, 
///     184.12117, 
///     271.1532, 
///     457.23251000000005,
///     186.07931, 
///     273.11134, 
///     344.14844999999997, 
///     457.23250999999993
/// ]);
/// assert_eq!(spectrum.unwrap(), result);
/// ```
pub fn complete_spectrum(seq: &Sequence) -> Result<Vec<f64>> {
    let mut spectrum = prefix_spectrum(seq)?;
    spectrum.extend(suffix_spectrum(seq)?);
    Ok(spectrum)
}

/// Infers the protein given from a given mass spectrum
///
/// # Arguments
/// * `seq` - sequence to generate complete spectrum
/// * `margin` - error tolerance when comparing spectrum masses with monoisotopic masses
///
/// # Example
/// ```
/// use biotech::analysis::spectro::infer_protein;
/// use biotech::ds::sequence::Sequence;
///
/// let spectrum = Vec::<f64>::from([
///     3524.8542, 3710.9335, 3841.974, 3970.0326, 4057.0646]);
/// let error_margin = 0.01_f64;
/// let inferred_protein = infer_protein(&spectrum, error_margin);
///
/// assert_eq!(inferred_protein.unwrap(), Sequence::from("WMQS"));
/// ```
pub fn infer_protein(spectrum: &[f64], margin: f64) -> Result<Sequence> {
    let mut result = Sequence::new();
    for i in 1..spectrum.len() {
        let diff = spectrum[i] - spectrum[i - 1];
        let x = monoisotopic_mass_table
            .iter()
            .find_map(|(key, &value)| {
                if (value - diff).abs() < margin {
                    Some(key)
                } else {
                    None
                }
            })
            .ok_or(BioError::ItemNotFound)?;
        result.push(*x);
    }
    Ok(result)
}

/// From a set of proteins, finds the one that minimizes the error between its spectrum and a reference spectrum.
///
/// # Arguments
/// * `tile` - container with all the proteins to evaluate
/// * `spectrum` - reference spectrum to compare with
///
/// # Example
/// ```
/// use biotech::analysis::spectro::match_protein;
/// use biotech::ds::sequence::Sequence;
/// use biotech::ds::tile::Tile;
///
/// let mut proteins = Tile::new();
/// 
/// proteins.extend(Vec::from([
///     Sequence::from("GSDMQS"),
///     Sequence::from("VWICN"), 
///     Sequence::from("IASWMQS"), 
///     Sequence::from("PVSMGAD")
/// ]));
///
/// let spectrum = Vec::<f64>::from([
///     445.17838, 
///     115.02694, 
///     186.07931, 
///     314.13789, 
///     317.1198, 
///     215.09061
/// ]);
///
/// let matched_protein = match_protein(&proteins, &spectrum);
/// assert_eq!(matched_protein.unwrap(), Sequence::from("PVSMGAD"));
/// ```
pub fn match_protein(tile: &Tile, spectrum: &[f64]) -> Result<Sequence> {
    let mut max = 0;
    let mut max_idx = 0;

    for (idx, s) in tile.into_iter().enumerate() {
        let ps = complete_spectrum(&s)?;
        let (cnt, _) = spectral_mass_shift(spectrum, &ps)?;
        if cnt >= max {
            max = cnt;
            max_idx = idx;
        }
    }

    Ok(tile[max_idx].clone())
}
