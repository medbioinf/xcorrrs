use context_error::BoxedError;
use mzannotate::fragment::Fragment;
use mzcore::{
    chemistry::{AmbiguousMolecule, MassMode, OutputMolecularFormula},
    ontology::Ontologies,
    sequence::PeptidoformIon,
};
use mzdata::prelude::SpectrumLike;
use mzpeaks::{CentroidPeak, DeconvolutedPeak};

use crate::{
    utils::{create_theoretical_fragments, extract_peaks, resolve_precursor_charge},
    xcorr::{configuration::FinalizedConfiguration, error::Error, scoring_result::ScoringResult},
};

/// +/- m/z shift for the xcorr calculation.
pub const BIN_SHIFT: usize = 75;

/// According to the original authors, after binning the experimental spectrum, the maximum intensity is normalized
/// over a number of fixed windows.
pub const NUM_WINDOWS_FOR_NORMALIZATION: u8 = 10;

/// Size of each chunk in the sparse SP score matrix.
const SP_MATRIX_SIZE: usize = 100;

/// A small value to consider as zero in the SP score matrix.
/// Called FLOAT_ZERO in the Comet source code.
const NEARLY_ZERO: f64 = 1e-6;

pub struct FastXcorr<'a> {
    config: &'a FinalizedConfiguration,
    fragment_charge: usize,
    /// Sparse SP score matrix
    sp_matrix: Vec<Option<Vec<f64>>>,
    /// y' prime from equation 6 in https://pubs.acs.org/doi/10.1021/pr800420s
    preprocessed_experimental_spectrum: Vec<f64>,
}

impl FastXcorr<'_> {
    /// Creates a new FastXcorr instance. This includes all preprocessing steps for the experimental spectrum.
    ///
    /// Arguments:
    /// * `config` - The configuration to use for scoring.
    /// * `experimental_spectrum` - The experimental spectrum to score against. Must report its
    ///   own precursor charge; scoring fails with [`Error::MissingCharge`] otherwise. Must also
    ///   report exactly one precursor with exactly one selected ion (e.g. not a chimeric
    ///   spectrum); scoring fails with [`Error::InvalidPrecursor`] otherwise.
    ///
    pub fn new<'a, S: SpectrumLike<CentroidPeak, DeconvolutedPeak>>(
        config: &'a FinalizedConfiguration,
        experimental_spectrum: &S,
    ) -> Result<FastXcorr<'a>, Error> {
        let charge = resolve_precursor_charge(experimental_spectrum)?.ok_or(Error::MissingCharge)?;

        let (raw_mz, raw_intensities) =
            extract_peaks(experimental_spectrum).map_err(|_| Error::UncentroidedSpectrum)?;
        let (mz, intensities): (Vec<f64>, Vec<f64>) = raw_mz
            .into_iter()
            .zip(raw_intensities)
            .filter(|&(mz, intensity)| {
                intensity >= config.minimum_intensity
                    && match config.clear_mz_range {
                        Some((min_mz, max_mz)) => mz <= min_mz && mz >= max_mz,
                        None => true,
                    }
            })
            .unzip();

        if mz.is_empty() {
            return Err(Error::EmptyExperimentalSpectrum);
        }

        // Get max m/z once for binning and normalization
        let mz_max = mz.iter().fold(f64::NEG_INFINITY, |acc, &x| acc.max(x));

        // Binning
        let binned_experimental_spectrum = Self::experimental_spectrum_binning(
            &mz,
            &intensities,
            mz_max,
            config.bin_size,
            config.bin_offset,
        )?;

        // Get max sqrt intensity for SP score matrix calculation
        let max_intensity_sqrt = binned_experimental_spectrum
            .iter()
            .fold(f64::NEG_INFINITY, |a, &b| a.max(b));

        let sp_matrix =
            Self::build_sparse_sp_score(&binned_experimental_spectrum, max_intensity_sqrt);

        // Normalization of the binned experimental spectrum
        let binned_normalized_experimental_spectrum = Self::experimental_spectrum_normalization(
            mz_max,
            &binned_experimental_spectrum,
            config.bin_size,
            config.bin_offset,
            config.use_flanking_peaks,
            max_intensity_sqrt,
        )?;

        let preprocessed_experimental_spectrum =
            Self::xcorr_prerprocessing(&binned_normalized_experimental_spectrum);

        let mut fragment_charge = (charge - 1).max(1);
        if fragment_charge > config.max_fragment_charge {
            fragment_charge = config.max_fragment_charge;
        }

        Ok(FastXcorr {
            config,
            fragment_charge,
            sp_matrix,
            preprocessed_experimental_spectrum,
        })
    }

    /// Calculates the number of bins based on the maximum m/z value, bin size, and bin offset.
    ///
    /// # Arguments:
    /// * `mz_max` - The maximum m/z value.
    /// * `bin_size` - The size of each bin.
    /// * `bin_offset` - The offset in monoisotopic mass.
    ///
    pub fn calc_number_of_bins(mz_max: f64, bin_size: f64, bin_offset: f64) -> usize {
        (mz_max / bin_size + 1.0 - bin_offset) as usize + 2 + BIN_SHIFT
    }

    /// Calculates the binned position for a given m/z value.
    /// This is used to determine which bin a given m/z value falls into.
    ///
    /// Arguments:
    /// * `mz` - The m/z value to be binned.
    /// * `bin_size` - The size of the bins to be used for binning the spectrum.
    /// * `bin_offset` - The offset in m/z to be applied before bin
    ///
    pub fn calc_binned_position(mz: f64, bin_size: f64, bin_offset: f64) -> usize {
        (mz / bin_size + 1.0 - bin_offset) as usize
    }

    /// Bins the experimental spectrum according to the provided bin size and offset.s
    ///
    /// # Arguments
    /// * `mz` - The m/z values of the experimental spectrum.
    /// * `intensities` - The intensities of the experimental spectrum.
    /// * `mz_max` - The maximum m/z value in the experimental spectrum.
    /// * `bin_size` - The size of the bins used for binning the experimental spectrum.
    /// * `bin_offset` - The offset in m/z to be applied before binning.
    ///
    pub fn experimental_spectrum_binning(
        mz: &[f64],
        intensities: &[f64],
        mz_max: f64,
        bin_size: f64,
        bin_offset: f64,
    ) -> Result<Vec<f64>, Error> {
        let number_of_bins = Self::calc_number_of_bins(mz_max, bin_size, bin_offset);
        let mut binned_spectrum: Vec<f64> = vec![0.0; number_of_bins];

        for (mz, intensity) in mz.iter().zip(intensities.iter()) {
            let index = Self::calc_binned_position(*mz, bin_size, bin_offset);
            binned_spectrum[index] = binned_spectrum[index].max(intensity.sqrt());
        }

        Ok(binned_spectrum)
    }

    /// Normalizes the binned experimental spectrum according to Comet's normalization procedure
    /// (`CometPreprocess::MakeCorrData`).
    ///
    /// # Arguments
    /// * `mz_max` - The maximum m/z value in the experimental spectrum.
    /// * `binned_experimental_spectrum` - The binned, raw (unnormalized) experimental spectrum.
    /// * `bin_size` - The size of the bins used for binning the experimental spectrum.
    /// * `bin_offset` - The offset in m/z to be applied before binning.
    /// * `use_flanking_peaks` - Whether to use flanking peaks (half the intensity on both sides) in the normalization.
    /// * `highest_intensity` - The maximum (sqrt) intensity across the whole binned spectrum, used
    ///   as the basis for the 5% noise cutoff (Comet compares against the spectrum-wide maximum,
    ///   not the local window maximum).
    ///
    pub fn experimental_spectrum_normalization(
        mz_max: f64,
        binned_experimental_spectrum: &[f64],
        bin_size: f64,
        bin_offset: f64,
        use_flanking_peaks: bool,
        highest_intensity: f64,
    ) -> Result<Vec<f64>, Error> {
        let highest_ion = Self::calc_binned_position(mz_max, bin_size, bin_offset);
        let windows_size = (highest_ion as f64 / NUM_WINDOWS_FOR_NORMALIZATION as f64) as usize + 1;
        // Comet only keeps bins that exceed 5% of the *spectrum-wide* maximum intensity; every
        // other bin is dropped to zero rather than left at its raw (unnormalized) value.
        let intensity_cutoff = 0.05 * highest_intensity;

        let mut normalized_spectrum = vec![0.0; binned_experimental_spectrum.len()];

        for window_start in (0..binned_experimental_spectrum.len()).step_by(windows_size) {
            let window_end = (window_start + windows_size).min(binned_experimental_spectrum.len());
            let window = &binned_experimental_spectrum[window_start..window_end];

            let window_max = window.iter().fold(f64::NEG_INFINITY, |acc, &x| acc.max(x));

            if window_max > 0.0 {
                let scale = 50.0 / window_max;
                for (offset, &value) in window.iter().enumerate() {
                    if value > intensity_cutoff {
                        normalized_spectrum[window_start + offset] = value * scale;
                    }
                }
            }
        }

        if !use_flanking_peaks {
            Ok(normalized_spectrum)
        } else {
            let mut flanked_binned_spectrum = vec![0.0; normalized_spectrum.len()];
            normalized_spectrum
                .into_iter()
                .enumerate()
                .for_each(|(i, value)| {
                    flanked_binned_spectrum[i] += value;
                    let half_peak = value * 0.5;

                    if i > 0 {
                        flanked_binned_spectrum[i - 1] += half_peak;
                    }
                    if i < flanked_binned_spectrum.len() - 1 {
                        flanked_binned_spectrum[i + 1] += half_peak;
                    }
                });
            Ok(flanked_binned_spectrum)
        }
    }

    /// Calculates y' from equation 6 in https://pubs.acs.org/doi/10.1021/pr800420s
    ///
    /// Arguments:
    /// * `binned_normalized_experimental_spectrum` - The binned and normalized experimental.
    ///
    pub fn xcorr_prerprocessing(binned_normalized_experimental_spectrum: &[f64]) -> Vec<f64> {
        // Extend by `shift` bins on both sides
        let mut corrected_experimental_spectrum_shift =
            vec![0.0; binned_normalized_experimental_spectrum.len()];

        // Shift -75 to -1  would be just a sum of the first 75 bins at position 0
        let mut sum_offsets: f64 = binned_normalized_experimental_spectrum[1..=BIN_SHIFT]
            .iter()
            .sum();
        let mean_offset = sum_offsets / 150.0;
        corrected_experimental_spectrum_shift[0] =
            binned_normalized_experimental_spectrum[0] - mean_offset;

        // Calculate this once to use it in the loop and make it depending on the BIN_SHIFT
        let bin_shift_plus = BIN_SHIFT + 1;

        // For each subsequent i, update the sliding window
        for i in 1..binned_normalized_experimental_spectrum.len() {
            if i >= bin_shift_plus {
                sum_offsets -= binned_normalized_experimental_spectrum[i - bin_shift_plus];
            }

            let add_idx = i + BIN_SHIFT;
            if add_idx < binned_normalized_experimental_spectrum.len() {
                sum_offsets += binned_normalized_experimental_spectrum[add_idx];
            }

            let old_center = i - 1;
            if old_center < binned_normalized_experimental_spectrum.len() {
                sum_offsets += binned_normalized_experimental_spectrum[old_center];
            }

            if i < binned_normalized_experimental_spectrum.len() {
                sum_offsets -= binned_normalized_experimental_spectrum[i];
            }

            let mean_offset = sum_offsets / 150.0;
            corrected_experimental_spectrum_shift[i] =
                binned_normalized_experimental_spectrum[i] - mean_offset;
        }

        corrected_experimental_spectrum_shift
    }

    /// Calculates the xcorr between an already binned theoretical and binned experimental spectra according to
    /// equation 6 in https://pubs.acs.org/doi/10.1021/pr800420s
    ///
    /// # Arguments
    /// * `theoretical_spectrum` - Theoretical fragments.
    /// * `preprocessed_experimental_spectrum` - The y' values calculated from the experimental spectrum.
    /// * `bin_size` - The size of the bins used for binning the experimental spectrum.
    /// * `bin_offset` - The offset in m/z to be applied before binning
    ///
    pub fn xcorr_spectra(
        theoretical_spectrum: &[f64],
        preprocessed_experimental_spectrum: &[f64],
        bin_size: f64,
        bin_offset: f64,
    ) -> f64 {
        let xcorr: f64 = theoretical_spectrum
            .iter()
            .map(|mz| {
                let index = Self::calc_binned_position(*mz, bin_size, bin_offset);
                if index < preprocessed_experimental_spectrum.len() {
                    preprocessed_experimental_spectrum[index]
                } else {
                    0.0
                }
            })
            .sum();

        xcorr * 0.005
    }

    /// Calculates the number of matched ions based on the sparse SP score matrix.
    ///
    /// # Arguments
    /// * `theoretical_spectrum` - Theoretical m/z
    /// * `sp_matrix` - Sparse SP score matrix
    /// * `bin_size` - Bin size for the spectra
    /// * `bin_offset` - Bin offset for the spectra
    ///
    fn match_ions_sp_score_based(
        theoretical_spectrum: &[f64],
        sp_matrix: &[Option<Vec<f64>>],
        bin_size: f64,
        bin_offset: f64,
    ) -> usize {
        theoretical_spectrum
            .iter()
            .map(|mz| {
                let index = Self::calc_binned_position(*mz, bin_size, bin_offset);
                let sp_score = Self::find_sp_score(sp_matrix, index);
                if sp_score > NEARLY_ZERO {
                    1
                } else {
                    0
                }
            })
            .sum()
    }

    /// Build a chunked sparse SP score matrix from a dense vector of intensities.
    ///
    /// Behavior mirrors the C implementation in `CometPreprocess::PreprocessSpectrum`:
    /// - Normalize intensities to `100.0 * intensity / max_intensity`.
    /// - SPARSE rows are `sparse_matrix_size` wide.
    /// - Number of rows `i_sp_score_data = data.len() / sparse_matrix_size + 1`.
    /// - Rows are allocated on-demand (only when at least one normalized value in the row > NEARLY_ZERO).
    /// - Returns `(sparse_matrix, i_sp_score_data)`, where `sparse_matrix` is a Vec<Option<Vec<f32>>>.
    ///
    /// # Arguments
    /// * `binned_theoretical_spectrum` - Binned theoretical spectrum.
    /// * `max_intensity_sqrt` - Maximum intensity in the binned normalized spectrum.
    ///
    pub fn build_sparse_sp_score(
        binned_theoretical_spectrum: &[f64],
        max_intensity_sqrt: f64,
    ) -> Vec<Option<Vec<f64>>> {
        // compute number of sparse rows (iSpScoreData)
        let matrix_size = binned_theoretical_spectrum.len() / SP_MATRIX_SIZE + 1;

        // top-level vector of optional rows, initially None (like new float*[iSpScoreData]())
        let mut sparse: Vec<Option<Vec<f64>>> = vec![None; matrix_size];

        // If max_intensity_sqrt <= 0.0, treat as all zeros -> nothing to allocate
        if max_intensity_sqrt <= 0.0 {
            return sparse;
        }

        for (i, peak) in binned_theoretical_spectrum.iter().enumerate() {
            let normalized = 100.0 * peak / max_intensity_sqrt;
            // treat small values as zero (C checks pfSpScoreData[i] > NEARLY_ZERO)
            if normalized > NEARLY_ZERO {
                let x = i / SP_MATRIX_SIZE; // row
                let y = i - (x * SP_MATRIX_SIZE); // col in row

                // allocate row if necessary, initialize to zeros
                if sparse[x].is_none() {
                    sparse[x] = Some(vec![0.0_f64; SP_MATRIX_SIZE]);
                }

                // safe to unwrap now
                if let Some(ref mut row) = sparse[x] {
                    row[y] = normalized;
                }
            }
        }

        sparse
    }

    /// Return the SP score for a given binned position from the sparse SP score matrix.
    ///
    /// # Arguments
    /// * `sparse` - The sparse SP score matrix.
    /// * `bin` - The binned position to retrieve the SP score for.
    ///
    pub fn find_sp_score(sparse: &[Option<Vec<f64>>], bin: usize) -> f64 {
        let x = bin / SP_MATRIX_SIZE; // row

        if x >= sparse.len() || bin == 0 || sparse[x].is_none() {
            return 0.0_f64;
        }

        let y = bin - (x * SP_MATRIX_SIZE); // col in row
        let row = sparse[x].as_ref().unwrap();
        row[y]
    }

    /// Calculates the theoretical fragments for a given peptide.
    ///
    pub fn create_theoretical_fragments(
        &self,
        peptidoform_ion: &PeptidoformIon,
    ) -> Result<Vec<Fragment<OutputMolecularFormula>>, Error> {
        create_theoretical_fragments(
            peptidoform_ion,
            &self.config.fragmentation_model,
            self.fragment_charge,
        )
        .map_err(|err| Error::TheoreticalSpectrum(Box::new(err)))
    }

    /// Creates an array of theoretical m/z values from a list of fragments.
    ///
    /// # Arguments
    /// * `theoretical_fragments` - Theoretical fragments to create the spectrum from.
    ///
    pub fn create_threoretical_mz(
        theoretical_fragments: &[Fragment<OutputMolecularFormula>],
    ) -> Vec<f64> {
        theoretical_fragments
            .iter()
            .filter_map(|f| f.mz(MassMode::Monoisotopic).map(|mz| mz.value))
            .collect::<Vec<f64>>()
    }

    /// Calculates the xcorr between a peptide and the experimental spectrum
    ///
    /// # Arguments
    /// * `peptide` - The peptide sequence to score.
    ///
    pub fn xcorr_peptide(
        &self,
        peptide: &str,
        ontologies: Option<&Ontologies>,
    ) -> Result<ScoringResult, Error> {
        let empty_ontologies = Ontologies::empty();
        let ontologies = ontologies.unwrap_or(&empty_ontologies);
        let (peptidoform_ion, parsing_errors) = PeptidoformIon::pro_forma(peptide, ontologies)
            .map_err(|err| {
                Error::InvalidProForma(err.into_iter().map(BoxedError::to_owned).collect())
            })?;

        if !parsing_errors.is_empty() {
            return Err(Error::InvalidProForma(
                parsing_errors
                    .into_iter()
                    .map(BoxedError::to_owned)
                    .collect(),
            ));
        }

        let (min_theoretical_mass, max_theoretical_mass) =
            match peptidoform_ion.formulas().mass_bounds().into_option() {
                Some((min, max)) => (min.monoisotopic_mass().value, max.monoisotopic_mass().value),
                None => (-1.0, -1.0),
            };

        let theoretical_fragments = self.create_theoretical_fragments(&peptidoform_ion)?;
        let theoretical_mz = Self::create_threoretical_mz(&theoretical_fragments);
        let ions_total = theoretical_mz.len();

        let ions_matched = Self::match_ions_sp_score_based(
            &theoretical_mz,
            &self.sp_matrix,
            self.config.bin_size,
            self.config.bin_offset,
        );

        // Score
        let score = Self::xcorr_spectra(
            &theoretical_mz,
            &self.preprocessed_experimental_spectrum,
            self.config.bin_size,
            self.config.bin_offset,
        );

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
    use std::{env, io::Write, path::PathBuf};

    use itertools::multiunzip;
    use polars::prelude::*;
    use rayon::prelude::*;

    use crate::{
        utils::tests::{build_raw_spectrum, get_spectrum, mean_sq_err},
        xcorr::{
            configuration::{Configuration, FinalizedConfiguration},
            fast_xcorr::FastXcorr,
        },
    };

    use super::*;

    /// Reads the test data from a TSV file and processes it to include a new column for proforma peptides.
    /// TEST_NUMBER_OF_PSMS environment variable limits to the X best scored PSMs.
    ///
    fn read_test_data() -> DataFrame {
        // Read the Comet results from a TSV file
        let mut comet_df = CsvReadOptions::default()
            .with_has_header(true)
            .with_parse_options(
                CsvParseOptions::default()
                    .with_separator(b'\t')
                    .with_comment_prefix(Some("#")),
            )
            .try_into_reader_with_file_path(Some(PathBuf::from(
                "test_files/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.comet.tsv",
            )))
            .unwrap()
            .finish()
            .unwrap();

        // Sort by xcorr in descending order
        comet_df
            .sort_in_place(
                ["xcorr"],
                SortMultipleOptions::default().with_order_descending(true),
            )
            .unwrap();

        // Reduce the number of PSMs if the environment variable is set
        let mut comet_df = match std::env::var("TEST_NUMBER_OF_PSMS") {
            Ok(number_of_psms) => comet_df.slice(0, number_of_psms.parse::<usize>().unwrap()),
            Err(_) => comet_df,
        };

        // Create the proforma peptides column
        let modified_peptide = comet_df.column("modified_peptide").unwrap().str().unwrap();
        let profoma_peptides = modified_peptide
            .iter()
            .filter_map(|s| {
                s.map(|s| {
                    let mut proform_string = s[2..s.len() - 2]
                        .to_string()
                        .replace("[15.9949]", "[+15.9949]");
                    if proform_string.contains("C") {
                        proform_string = proform_string.replace("C", "C[+57.021464]");
                    }
                    proform_string
                })
            })
            .collect::<Vec<String>>();

        comet_df
            .with_column(
                Series::new("proforma_peptide".into(), profoma_peptides)
                    .cast(&DataType::String)
                    .unwrap(),
            )
            .unwrap();

        comet_df
    }

    /// Reads the experimental spectrum Eng provided for the peptide DIGSETK.
    ///
    pub fn get_eng_experimental_spectrum() -> (Vec<f64>, Vec<f64>) {
        let spec_df =
            ParquetReader::new(std::fs::File::open("test_files/eng/DIGSETK.parquet").unwrap())
                .read_parallel(ParallelStrategy::None)
                .finish()
                .unwrap();

        let mz_array = spec_df["mz"]
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<f64>>();

        let intensity_array = spec_df["intensity"]
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<f64>>();

        (mz_array, intensity_array)
    }

    /// Reads the partial xcorr spectrum and the index provided by Eng.
    /// Return is a tuple with the index and the fast xcorr values.
    ///
    pub fn get_eng_fast_xcorr_spectrum() -> (Vec<usize>, Vec<f64>) {
        let spec_df = CsvReadOptions::default()
            .with_has_header(true)
            .with_rechunk(true)
            .with_parse_options(
                CsvParseOptions::default()
                    .with_separator(b'\t')
                    .with_comment_prefix(Some("#")),
            )
            .try_into_reader_with_file_path(Some(PathBuf::from(
                "test_files/eng/DIGSETK.process.tsv",
            )))
            .unwrap()
            .finish()
            .unwrap();

        (
            spec_df
                .column("index")
                .unwrap()
                .i64()
                .unwrap()
                .into_no_null_iter()
                .map(|x| x as usize)
                .collect::<Vec<usize>>(),
            spec_df
                .column("fast_xcorr")
                .unwrap()
                .f64()
                .unwrap()
                .into_no_null_iter()
                .collect::<Vec<f64>>(),
        )
    }

    /// A spectrum with more than one selected ion on its precursor (e.g. a chimeric spectrum)
    /// must be rejected, not silently scored against the first ion.
    #[test]
    fn test_new_rejects_multiple_selected_ions() {
        let config: FinalizedConfiguration = Configuration::default().into();
        let mut spectrum = build_raw_spectrum(&[100.0], &[1.0], 2, 300.0);
        spectrum.description.precursor[0]
            .ions
            .push(mzdata::spectrum::SelectedIon {
                mz: 305.0,
                charge: Some(3),
                ..Default::default()
            });

        assert!(matches!(
            FastXcorr::new(&config, &spectrum),
            Err(Error::InvalidPrecursor(
                crate::utils::PrecursorMultiplicityError::MultipleSelectedIons(2)
            ))
        ));
    }

    /// Checks the Xcorr spectrum (y-prime) against the data provided by Eng
    ///
    #[test]
    fn test_preprocessed_experimental_spectrum() {
        let expected_xcorr_spec = get_eng_fast_xcorr_spectrum();

        // Create config  for low resolution data
        let config: FinalizedConfiguration = Configuration {
            bin_size: 1.0005,
            bin_offset: 0.4,
            use_flanking_peaks: true,
            ..Configuration::default()
        }
        .into();

        let experimental_spectrum = get_eng_experimental_spectrum();
        let spectrum =
            build_raw_spectrum(&experimental_spectrum.0, &experimental_spectrum.1, 1, 0.0);

        let xcorr = FastXcorr::new(&config, &spectrum).unwrap();

        let rouneded_xcorr_sped = xcorr
            .preprocessed_experimental_spectrum
            .iter()
            .map(|x| (x * 100.0).round() / 100.0)
            .collect::<Vec<f64>>();

        // Just select the values that are in the expected spectrum from Eng
        let rouneded_xcorr_sped = expected_xcorr_spec
            .0
            .iter()
            .map(|&index| rouneded_xcorr_sped[index])
            .collect::<Vec<f64>>();

        let rounded_expected_xcorr_spec = expected_xcorr_spec
            .1
            .iter()
            .map(|x| (x * 100.0).round() / 100.0)
            .collect::<Vec<f64>>();

        if env::var("VERBOSE").is_ok() {
            let max_mz = experimental_spectrum
                .0
                .iter()
                .fold(f64::NEG_INFINITY, |acc, &x| acc.max(x));

            let binned_experimental_spectrum = FastXcorr::experimental_spectrum_binning(
                &experimental_spectrum.0,
                &experimental_spectrum.1,
                max_mz,
                config.bin_size,
                config.bin_offset,
            )
            .unwrap();

            let highest_intensity = binned_experimental_spectrum
                .iter()
                .fold(f64::NEG_INFINITY, |a, &b| a.max(b));

            let binned_normalized_experimental_spectrum =
                FastXcorr::experimental_spectrum_normalization(
                    max_mz,
                    &binned_experimental_spectrum,
                    config.bin_size,
                    config.bin_offset,
                    config.use_flanking_peaks,
                    highest_intensity,
                )
                .unwrap();

            let binned_normalized_experimental_spectrum = expected_xcorr_spec
                .0
                .iter()
                .map(|&index| binned_normalized_experimental_spectrum[index])
                .collect::<Vec<f64>>();

            let output_file =
                std::fs::File::create("DIGSETK_fast_xcorr_preprocessing.tsv").unwrap();
            let mut output_writer = std::io::BufWriter::new(output_file);

            let _ = output_writer
                .write("bin\texpected\tcalc_binned\tcalc_prep\n".as_bytes())
                .unwrap();
            for (idx, (expected, (calc_binned, calc_prep))) in rounded_expected_xcorr_spec
                .iter()
                .zip(
                    binned_normalized_experimental_spectrum
                        .iter()
                        .zip(rouneded_xcorr_sped.iter()),
                )
                .enumerate()
            {
                let _ = output_writer
                    .write(format!("{idx}\t{expected}\t{calc_binned}\t{calc_prep}\n").as_bytes())
                    .unwrap();
            }
        }

        assert_eq!(rouneded_xcorr_sped, rounded_expected_xcorr_spec)
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
                .into_no_null_iter()
                .collect::<Vec<f64>>(),
            experimental_spectrum["intensity"]
                .f64()
                .unwrap()
                .into_no_null_iter()
                .collect::<Vec<f64>>(),
        );

        let config: FinalizedConfiguration = Configuration {
            bin_size: 1.0005,
            bin_offset: 0.4,
            ..Configuration::default()
        }
        .into();

        let spectrum =
            build_raw_spectrum(&experimental_spectrum.0, &experimental_spectrum.1, 1, 0.0);

        let xcorr = FastXcorr::new(&config, &spectrum).unwrap();

        if env::var("VERBOSE").is_ok() {
            let peptide = PeptidoformIon::pro_forma("DIGSETK", &Ontologies::empty())
                .unwrap()
                .0;

            let mut binned_theoretical_spectrum =
                vec![0.0; xcorr.preprocessed_experimental_spectrum.len()];

            for fragment in &xcorr.create_theoretical_fragments(&peptide).unwrap() {
                let mz = fragment.mz(MassMode::Monoisotopic).unwrap().value;
                let index = FastXcorr::calc_binned_position(mz, config.bin_size, config.bin_offset);
                if index < binned_theoretical_spectrum.len() {
                    binned_theoretical_spectrum[index] = 50.0;
                }
            }

            let output_file =
                std::fs::File::create("DIGSETK__fast_xcorr_bin___theoretical_bin.tsv").unwrap();
            let mut output_writer = std::io::BufWriter::new(output_file);

            let _ = output_writer
                .write("bin\texperimental_bin\ttheoretical_bin\n".as_bytes())
                .unwrap();
            for (idx, (fast_xcorr_bin, theoretical_bin)) in xcorr
                .preprocessed_experimental_spectrum
                .iter()
                .zip(binned_theoretical_spectrum.iter())
                .enumerate()
            {
                let _ = output_writer
                    .write(format!("{idx}\t{fast_xcorr_bin}\t{theoretical_bin}\n").as_bytes())
                    .unwrap();
            }
        }

        let scoring = xcorr.xcorr_peptide("DIGSETK", None).unwrap();
        println!("{scoring}");
        assert_eq!((scoring.score * 100.0).round() / 100.0, 2.92);
    }

    #[test]
    fn test_xcorr() {
        let comet_df = read_test_data();

        #[allow(clippy::type_complexity)]
        let results: Vec<(i64, String, f64, f64, u64, u64, u64, u64)> = (0..comet_df.height())
            .into_par_iter()
            .map(|idx| {
                let scan = comet_df["scan"].i64().unwrap().get(idx).unwrap();
                let comet_xcorr = comet_df["xcorr"].f64().unwrap().get(idx).unwrap();
                let comet_ions_total =
                    comet_df["ions_total"].i64().unwrap().get(idx).unwrap() as u64;
                let comet_ions_matched =
                    comet_df["ions_matched"].i64().unwrap().get(idx).unwrap() as u64;
                let proforma_peptide = comet_df["proforma_peptide"]
                    .str()
                    .unwrap()
                    .get(idx)
                    .unwrap();

                let spectrum = get_spectrum(scan.to_string().as_str());

                let config: FinalizedConfiguration = Configuration {
                    use_flanking_peaks: true,
                    common: crate::utils::FragmentIonConfiguration {
                        max_fragment_charge: 5,
                        ..crate::utils::FragmentIonConfiguration::default()
                    },
                    ..Configuration::default()
                }
                .into();

                // fast xcorr implementation
                let fast_xcorr = FastXcorr::new(&config, &spectrum).unwrap();

                let scoring = fast_xcorr.xcorr_peptide(proforma_peptide, None).unwrap();

                (
                    scan,
                    proforma_peptide.to_string(),
                    comet_xcorr,
                    scoring.round_score(3),
                    comet_ions_matched,
                    scoring.ions_matched as u64,
                    comet_ions_total,
                    scoring.ions_total as u64,
                )
            })
            .collect();

        let (
            scan_col,
            peptide_col,
            comet_xcorr_col,
            rs_xcorr_col,
            comet_ions_match_col,
            rs_ions_matched_col,
            comet_ions_total_col,
            rs_ions_total_col,
        ): (
            Vec<_>,
            Vec<_>,
            Vec<_>,
            Vec<_>,
            Vec<_>,
            Vec<_>,
            Vec<_>,
            Vec<_>,
        ) = multiunzip(results);

        let mut xcorrrs_df = DataFrame::new(vec![
            Column::new("scan".into(), scan_col),
            Column::new("modified_peptide".into(), peptide_col),
            Column::new("comet_xcorr".into(), comet_xcorr_col),
            Column::new("rs_xcorr".into(), rs_xcorr_col),
            Column::new("comet_ions_match".into(), comet_ions_match_col),
            Column::new("rs_ions_matched".into(), rs_ions_matched_col),
            Column::new("comet_ions_total".into(), comet_ions_total_col),
            Column::new("rs_ions_total".into(), rs_ions_total_col),
        ])
        .unwrap();

        CsvWriter::new(std::fs::File::create("comparison.tsv").unwrap())
            .with_separator(b'\t')
            .finish(&mut xcorrrs_df)
            .unwrap();

        let comet_xcorr_max = xcorrrs_df
            .column("comet_xcorr")
            .unwrap()
            .f64()
            .unwrap()
            .max()
            .unwrap();

        let fast_xcorrrs_max = xcorrrs_df
            .column("rs_xcorr")
            .unwrap()
            .f64()
            .unwrap()
            .max()
            .unwrap();

        let max_score = comet_xcorr_max.max(fast_xcorrrs_max);

        if env::var("VERBOSE").is_ok() {
            let mut plot = plotly::Plot::new();
            let diagonal_trace =
                plotly::Scatter::new(vec![0.0, comet_xcorr_max], vec![0.0, fast_xcorrrs_max])
                    .mode(plotly::common::Mode::Lines)
                    .marker(plotly::common::Marker::default().color("red"))
                    .hover_info(plotly::common::HoverInfo::None)
                    .show_legend(false);

            let correlation_trace = plotly::Scatter::new(
                xcorrrs_df["comet_xcorr"].f64().unwrap().to_vec(),
                xcorrrs_df["rs_xcorr"].f64().unwrap().to_vec(),
            )
            .hover_text_array(
                xcorrrs_df["scan"]
                    .i64()
                    .unwrap()
                    .iter()
                    .map(|scan_opt| {
                        scan_opt
                            .map(|scan| format!("scan = {}", scan))
                            .unwrap_or_else(|| "scan = N/A".to_string())
                    })
                    .collect::<Vec<String>>(),
            )
            .mode(plotly::common::Mode::Markers)
            .marker(plotly::common::Marker::default().color("blue"))
            .show_legend(false);

            plot.add_trace(diagonal_trace);
            plot.add_trace(correlation_trace);

            plot.set_layout(
                plotly::Layout::new()
                    .title("Comet xcorr vs mzscoring xcorr")
                    .x_axis(
                        plotly::layout::Axis::new()
                            .title("Comet xcorr")
                            .constrain(plotly::layout::AxisConstrain::Domain),
                    )
                    .y_axis(
                        plotly::layout::Axis::new()
                            .title("mzscoring xcorr")
                            .scale_anchor("x"),
                    ),
            );
            plot.write_html("99-fast_xcorrrs_vs_comet_xcorr.html");
        }

        // Normalize comet xcorrs and calculates xcorrs

        let scaled_comet_xcorr = xcorrrs_df
            .column("comet_xcorr")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .map(|x| x / max_score)
            .collect::<Vec<f64>>();

        let scaled_fast_xcorrrs = xcorrrs_df
            .column("rs_xcorr")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .map(|x| x / max_score)
            .collect::<Vec<f64>>();

        let rmse_fast_xcorrs = mean_sq_err(&scaled_comet_xcorr, &scaled_fast_xcorrrs);

        println!("RMSE comet xcorr vs fast xcorrrs: {rmse_fast_xcorrs}");

        assert!(
            rmse_fast_xcorrs < 0.0002,
            "fast xcorr RMSE {rmse_fast_xcorrs} >= 0.0002"
        );
    }
}
