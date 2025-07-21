use ndarray::{s, Array1, Axis};
use rustyms::CompoundPeptidoformIon;

use crate::{
    binning::{experimental_spectrum_binning, theoretical_spectrum_binning},
    configuration::Configuration,
    error::Error,
    scoring_result::ScoringResult,
    utils::{calculate_number_of_bins_to_shift, create_threoretical_spectrum},
};

pub struct FastXcorr<'a> {
    config: &'a Configuration,
    charge: usize,
    filtered_experimental_spectrum: (Array1<f64>, Array1<f64>),
    max_experimental_mz: f64,
    fragment_charge: usize,
    shift: usize,
    y_prime: Array1<f64>,
}

impl FastXcorr<'_> {
    /// Creates a new FastXcorr instance.
    ///
    /// Arguments:
    /// * `config` - The configuration to use for scoring.
    /// * `experimental_spectrum` - The experimental spectrum to score against.
    /// * `charge` - Precursor charge
    ///
    pub fn new<'a>(
        config: &'a Configuration,
        experimental_spectrum: (&'a Array1<f64>, &'a Array1<f64>),
        charge: usize,
    ) -> Result<FastXcorr<'a>, Error> {
        if experimental_spectrum.0.is_empty() {
            return Err(Error::EmptyExperimentalSpectrum);
        }

        if experimental_spectrum.0.len() != experimental_spectrum.1.len() {
            return Err(Error::ExperimentalSpectrumShape(
                experimental_spectrum.0.len(),
                experimental_spectrum.1.len(),
            ));
        }

        // Filter out peaks below the minimum intensity
        let considerable_peaks_indexes = experimental_spectrum
            .1
            .iter()
            .enumerate()
            .filter(|(_, &intensity)| intensity >= config.minimum_intensity)
            .map(|(index, _)| index)
            .collect::<Vec<usize>>();

        let filtered_experimental_spectrum = (
            experimental_spectrum
                .0
                .select(Axis(0), &considerable_peaks_indexes),
            experimental_spectrum
                .1
                .select(Axis(0), &considerable_peaks_indexes),
        );

        let max_experimental_mz = *filtered_experimental_spectrum.0.last().unwrap();

        let shift = calculate_number_of_bins_to_shift(config.bin_size);
        let binned_experimental_spectrum = experimental_spectrum_binning(
            &filtered_experimental_spectrum.0,
            &filtered_experimental_spectrum.1,
            config.bin_size,
            config.bin_offset,
            charge,
            shift,
            config.use_flanking_peaks,
        )?;

        let y_prime = Self::calc_y_prime(&binned_experimental_spectrum, shift);

        let mut fragment_charge = (charge - 1).max(1);
        if fragment_charge > config.max_fragment_charge {
            fragment_charge = config.max_fragment_charge;
        }

        Ok(FastXcorr {
            config,
            charge,
            filtered_experimental_spectrum,
            max_experimental_mz,
            fragment_charge,
            shift,
            y_prime,
        })
    }

    /// Calculates the specturm +/- 75 m/z shift from equation 6 in https://pubs.acs.org/doi/10.1021/pr800420s
    ///
    /// Arguments:
    /// * `binned_experimental_spectrum` - The binned experimental m/z values +/- m/z shift.
    /// * `shift` - Number of shifted bins applied to the binned experimental spectrum.
    ///
    pub fn calc_y_prime_shift(
        binned_experimental_spectrum: &Array1<f64>,
        shift: usize,
    ) -> Array1<f64> {
        // Extend by `shift` bins on both sides
        let mut y_prime_shift = Array1::zeros(binned_experimental_spectrum.len());

        // Binned spectrum wihtout the shifts
        let shiftless_binned_spectrum = binned_experimental_spectrum
            .slice(s![shift..binned_experimental_spectrum.len() - shift]);

        // Shift -75 to -1
        for slice_start in 0..shift {
            let slice_end = slice_start + shiftless_binned_spectrum.len();
            let mut y_prime_shift_slice = y_prime_shift.slice_mut(s![slice_start..slice_end]);
            y_prime_shift_slice += &shiftless_binned_spectrum;
        }

        // Shift +1 to +75
        let shift_start = shift + 1;
        let shift_end = shift_start + shift;
        for slice_start in shift_start..shift_end {
            let slice_end = slice_start + shiftless_binned_spectrum.len();
            let mut y_prime_shift_slice = y_prime_shift.slice_mut(s![slice_start..slice_end]);
            y_prime_shift_slice += &shiftless_binned_spectrum;
        }

        y_prime_shift
    }

    /// Calculates y' from equation 6 in https://pubs.acs.org/doi/10.1021/pr800420s
    ///
    /// Arguments:
    /// * `binned_experimental_spectrum` - The binned experimental m/z values +/- m/z shift.
    /// * `shift` - Number of shifted bins applied to the binned experimental spectrum.
    ///
    pub fn calc_y_prime(binned_experimental_spectrum: &Array1<f64>, shift: usize) -> Array1<f64> {
        let mut y_prime_shift = Self::calc_y_prime_shift(binned_experimental_spectrum, shift);
        y_prime_shift /= 150.0;
        binned_experimental_spectrum - &y_prime_shift
    }

    /// Calculates the xcorr between an already binned theoretical and binned experimental spectra according to
    /// equation 6 in https://pubs.acs.org/doi/10.1021/pr800420s
    ///
    /// # Arguments
    /// * `binned_theoretical_spectrum` - The binned theoretical m/z values + 75 m/z shift.
    /// * `y_prime` - The y' values calculated from the experimental spectrum.
    ///
    pub fn xcorr_binned_spectrum(
        binned_theoretical_spectrum: &Array1<f64>,
        y_prime: &Array1<f64>,
    ) -> f64 {
        binned_theoretical_spectrum.dot(y_prime) / 10000.0
    }

    pub fn create_threoretical_spectrum(
        &self,
        peptide: &CompoundPeptidoformIon,
    ) -> Result<Array1<f64>, Error> {
        create_threoretical_spectrum(
            peptide,
            &self.config.fragmentation_model,
            self.fragment_charge,
            self.max_experimental_mz,
        )
    }

    pub fn theoretical_spectrum_binning(
        &self,
        theoretical_spectrum: &Array1<f64>,
    ) -> Result<Array1<f64>, Error> {
        theoretical_spectrum_binning(
            theoretical_spectrum,
            self.config.bin_size,
            self.config.bin_offset,
            self.charge,
            self.shift,
            Some(self.max_experimental_mz),
            false, // In the fast xcorr implementation the flanking peaks where "moved" to the experimental spectrum
        )
    }

    /// Calculates the xcorr between a peptide and the experimental spectrum
    ///
    /// # Arguments
    /// * `peptide` - The peptide sequence to score.
    ///
    pub fn xcorr_peptide(&self, peptide: &str) -> Result<ScoringResult, Error> {
        let peptide = CompoundPeptidoformIon::pro_forma(peptide, None)
            .map_err(Error::InvalidPeptideSequence)?;

        let (min_theoretical_mass, max_theoretical_mass) =
            match peptide.formulas().mass_bounds().into_option() {
                Some((min, max)) => (min.monoisotopic_mass().value, max.monoisotopic_mass().value),
                None => (-1.0, -1.0),
            };

        let theoretical_spectrum = self.create_threoretical_spectrum(&peptide)?;
        let ions_total = theoretical_spectrum.len();

        let bin_offset_mz = if self.config.bin_offset > 0.0 {
            crate::utils::dalton_to_mass_to_charge(self.config.bin_offset, self.charge)
        } else {
            0.0
        };
        let ions_matched = self
            .filtered_experimental_spectrum
            .0
            .iter()
            .zip(self.filtered_experimental_spectrum.1.iter())
            .map(|(&mz, &intensity)| {
                let index =
                    ((mz + bin_offset_mz) / self.config.bin_size).floor() as usize + self.shift - 1;
                if self.y_prime[index] > 0.0 && intensity > 0.0 {
                    1
                } else {
                    0
                }
            })
            .sum();

        let binned_thereoretical_spectrum =
            self.theoretical_spectrum_binning(&theoretical_spectrum)?;

        drop(theoretical_spectrum);

        // Score
        let score = Self::xcorr_binned_spectrum(&binned_thereoretical_spectrum, &self.y_prime);

        // binned_thereoretical_spectrum
        //     .slice_axis_inplace(Axis(0), Slice::new((self.shift + 40) as isize, None, 1));

        Ok(ScoringResult {
            score,
            min_theoretical_mass,
            max_theoretical_mass,
            ions_total,
            ions_matched,
        })
    }
}

#[cfg(test)]
mod tests {
    use std::{env, io::Write};

    use crate::utils::tests::{get_eng_experimental_spectrum, get_eng_fast_xcorr_spectrum};

    use ndarray::Slice;
    use polars::prelude::*;

    use super::*;

    /// Checks the Xcorr spectrum (y-prime) against the data provided by Eng
    ///
    #[test]
    fn test_y_prime() {
        let expected_xcorr_spec = get_eng_fast_xcorr_spectrum();

        // Create config  for low resolution data
        let config = Configuration {
            bin_size: 1.0005,
            bin_offset: 0.4,
            ..Configuration::default()
        };

        let experimental_spectrum = get_eng_experimental_spectrum();

        let xcorr = FastXcorr::new(
            &config,
            (&experimental_spectrum.0, &experimental_spectrum.1),
            1,
        )
        .unwrap();

        let mut rouneded_xcorr_sped = xcorr
            .y_prime
            .iter()
            .map(|x| (x * 100.0).round() / 100.0)
            .collect::<Array1<f64>>();

        // Remove the front shift
        rouneded_xcorr_sped.slice_axis_inplace(Axis(0), Slice::new(xcorr.shift as isize, None, 1));

        // Just select the values that are in the expected spectrum from Eng
        rouneded_xcorr_sped =
            rouneded_xcorr_sped.select(Axis(0), expected_xcorr_spec.0.as_slice().unwrap());

        let rounded_expeced_xcorr_spec = expected_xcorr_spec
            .1
            .iter()
            .map(|x| (x * 100.0).round() / 100.0)
            .collect::<Array1<f64>>();

        assert_eq!(rouneded_xcorr_sped, rounded_expeced_xcorr_spec)
    }

    /// Tests the xcorr calculation against data provided by the J. Eng
    ///
    #[test]
    fn test_xcorr_eng_data() {
        // Load experimental spectrum from Parquet file
        let experimental_spectrum =
            ParquetReader::new(std::fs::File::open("test_files/eng/DIGSETK.parquet").unwrap())
                .read_parallel(ParallelStrategy::None)
                .finish()
                .unwrap();

        let experimental_spectrum = (
            experimental_spectrum["mz"]
                .f64()
                .unwrap()
                .to_ndarray()
                .unwrap()
                .to_owned(),
            experimental_spectrum["intensity"]
                .f64()
                .unwrap()
                .to_ndarray()
                .unwrap()
                .to_owned(),
        );

        let config = Configuration {
            bin_size: 1.0005,
            bin_offset: 0.4,
            use_flanking_peaks: true,
            ..Configuration::default()
        };

        let xcorr = FastXcorr::new(
            &config,
            (&experimental_spectrum.0, &experimental_spectrum.1),
            1,
        )
        .unwrap();

        if env::var("VERBOSE").is_ok() {
            let peptide = CompoundPeptidoformIon::pro_forma("DIGSETK", None).unwrap();
            let binned_theoretical_spectrum = xcorr
                .theoretical_spectrum_binning(
                    &xcorr.create_threoretical_spectrum(&peptide).unwrap(),
                )
                .unwrap();

            let output_file =
                std::fs::File::create("DIGSETK__fast_xcorr_bin___theoretical_bin.tsv").unwrap();
            let mut output_writer = std::io::BufWriter::new(output_file);

            let _ = output_writer
                .write("bin\texperimental_bin\ttheoretical_bin\n".as_bytes())
                .unwrap();
            for (idx, (fast_xcorr_bin, theoretical_bin)) in xcorr
                .y_prime
                .iter()
                .zip(binned_theoretical_spectrum.iter())
                .enumerate()
            {
                let _ = output_writer
                    .write(format!("{idx}\t{fast_xcorr_bin}\t{theoretical_bin}\n").as_bytes())
                    .unwrap();
            }
        }

        let scoring = xcorr.xcorr_peptide("DIGSETK").unwrap();
        println!("{scoring}");
        assert_eq!((scoring.score * 100.0).round() / 100.0, 2.92);
    }
}
