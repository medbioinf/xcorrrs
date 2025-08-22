use ndarray::{s, Array1, Axis};
use rustyms::CompoundPeptidoformIon;

use crate::{
    configuration::FinalizedConfiguration, error::Error, scoring_result::ScoringResult,
    utils::create_threoretical_spectrum,
};

/// +/- m/z shift for the xcorr calculation.
pub const BIN_SHIFT: usize = 75;

/// According to the original authors, after binning the experimental spectrum, the maximum intensity is normalized
/// over a number of fixed windows.
pub const NUM_WINDOWS_FOR_NORMALIZATION: u8 = 10;

pub struct FastXcorr<'a> {
    config: &'a FinalizedConfiguration,
    max_experimental_mz: f64,
    fragment_charge: usize,
    /// y' prime from equation 6 in https://pubs.acs.org/doi/10.1021/pr800420s
    preprocessed_experimental_spectrum: Array1<f64>,
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
        config: &'a FinalizedConfiguration,
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

        let mut filtered_experimental_spectrum = (
            experimental_spectrum
                .0
                .select(Axis(0), &considerable_peaks_indexes),
            experimental_spectrum
                .1
                .select(Axis(0), &considerable_peaks_indexes),
        );

        if let Some((min_mz, max_mz)) = config.clear_mz_range {
            let considerable_peaks_indexes = experimental_spectrum
                .0
                .iter()
                .enumerate()
                .filter(|(_, &mz)| mz <= min_mz && mz >= max_mz)
                .map(|(index, _)| index)
                .collect::<Vec<usize>>();

            filtered_experimental_spectrum = (
                experimental_spectrum
                    .0
                    .select(Axis(0), &considerable_peaks_indexes),
                experimental_spectrum
                    .1
                    .select(Axis(0), &considerable_peaks_indexes),
            );
        }

        let max_experimental_mz = *filtered_experimental_spectrum.0.last().unwrap();

        let binned_experimental_spectrum = Self::experimental_spectrum_binning(
            &filtered_experimental_spectrum.0,
            &filtered_experimental_spectrum.1,
            config.bin_size,
            config.bin_offset,
            config.use_flanking_peaks,
        )?;

        let preprocessed_experimental_spectrum =
            Self::xcorr_prerprocessing(&binned_experimental_spectrum);

        let mut fragment_charge = (charge - 1).max(1);
        if fragment_charge > config.max_fragment_charge {
            fragment_charge = config.max_fragment_charge;
        }

        Ok(FastXcorr {
            config,
            max_experimental_mz,
            fragment_charge,
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

    pub fn experimental_spectrum_binning(
        mz: &Array1<f64>,
        intensities: &Array1<f64>,
        bin_size: f64,
        bin_offset: f64,
        use_flanking_peaks: bool,
    ) -> Result<Array1<f64>, Error> {
        let mz_max = mz.iter().fold(f64::NEG_INFINITY, |acc, &x| acc.max(x));
        let number_of_bins = Self::calc_number_of_bins(mz_max, bin_size, bin_offset);
        let mut binned_spectrum: Array1<f64> = Array1::zeros(number_of_bins);

        for (mz, intensity) in mz.iter().zip(intensities.iter()) {
            let index = Self::calc_binned_position(*mz, bin_size, bin_offset);
            binned_spectrum[index] = binned_spectrum[index].max(intensity.sqrt());
        }

        let highest_ion = Self::calc_binned_position(mz_max, bin_size, bin_offset);
        let windows_size = (highest_ion as f64 / NUM_WINDOWS_FOR_NORMALIZATION as f64) as usize + 1;

        for window_start in (0..binned_spectrum.len()).step_by(windows_size) {
            let window_end = (window_start + windows_size).min(binned_spectrum.len());
            let mut window = binned_spectrum.slice_mut(s![window_start..window_end]);

            let window_max = window.iter().fold(f64::NEG_INFINITY, |acc, &x| acc.max(x));
            // Comet is only normalizing if the value is greater than 0.05 of the window's maximum intensity
            let window_intensity_cutoff = 0.05 * window_max;
            window.mapv_inplace(|value| {
                if value <= window_intensity_cutoff {
                    value
                } else {
                    value / window_max * 50.0
                }
            });
        }

        if !use_flanking_peaks {
            Ok(binned_spectrum)
        } else {
            let mut flanked_binned_spectrum = Array1::zeros(binned_spectrum.len());
            binned_spectrum
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
    /// * `binned_experimental_spectrum` - The binned experimental m/z values + m/z shift.
    ///
    pub fn xcorr_prerprocessing(binned_experimental_spectrum: &Array1<f64>) -> Array1<f64> {
        // Extend by `shift` bins on both sides
        let mut corrected_experimental_spectrum_shift =
            Array1::zeros(binned_experimental_spectrum.len());

        // Shift -75 to -1  would be just a sum of the first 75 bins at position 0
        let mut sum_offsets = binned_experimental_spectrum.slice(s![1..=BIN_SHIFT]).sum();
        let mean_offset = sum_offsets / 150.0;
        corrected_experimental_spectrum_shift[0] = binned_experimental_spectrum[0] - mean_offset;

        // Calculate this once to use it in the loop and make it depending on the BIN_SHIFT
        let bin_shift_plus = BIN_SHIFT + 1;

        // For each subsequent i, update the sliding window
        for i in 1..binned_experimental_spectrum.len() {
            if i >= bin_shift_plus {
                sum_offsets -= binned_experimental_spectrum[i - bin_shift_plus];
            }

            let add_idx = i + BIN_SHIFT;
            if add_idx < binned_experimental_spectrum.len() {
                sum_offsets += binned_experimental_spectrum[add_idx];
            }

            let old_center = i - 1;
            if old_center < binned_experimental_spectrum.len() {
                sum_offsets += binned_experimental_spectrum[old_center];
            }

            if i < binned_experimental_spectrum.len() {
                sum_offsets -= binned_experimental_spectrum[i];
            }

            let mean_offset = sum_offsets / 150.0;
            corrected_experimental_spectrum_shift[i] =
                binned_experimental_spectrum[i] - mean_offset;
        }

        corrected_experimental_spectrum_shift
    }

    /// Calculates the xcorr between an already binned theoretical and binned experimental spectra according to
    /// equation 6 in https://pubs.acs.org/doi/10.1021/pr800420s
    ///
    /// # Arguments
    /// * `theoretical_spectrum` - Theoretical fragments.
    /// * `preprocessed_experimental_spectrum` - The y' values calculated from the experimental spectrum.
    ///
    pub fn xcorr_spectra(
        theoretical_spectrum: &Array1<f64>,
        preprocessed_experimental_spectrum: &Array1<f64>,
        bin_size: f64,
        bin_offset: f64,
    ) -> f64 {
        let xcorr: f64 = theoretical_spectrum
            .iter()
            .map(|mz| {
                let index = Self::calc_binned_position(*mz, bin_size, bin_offset);
                preprocessed_experimental_spectrum[index]
            })
            .sum();

        xcorr * 0.005
    }

    /// Calculates the xcorr between an already binned theoretical and binned experimental spectra according to
    /// equation 6 in https://pubs.acs.org/doi/10.1021/pr800420s
    ///
    /// # Arguments
    /// * `theoretical_spectrum` - Theoretical fragments.
    /// * `preprocessed_experimental_spectrum` - The y' values calculated from the experimental spectrum.
    ///
    pub fn matched_ions(
        theoretical_spectrum: &Array1<f64>,
        preprocessed_experimental_spectrum: &Array1<f64>,
        bin_size: f64,
        bin_offset: f64,
    ) -> usize {
        theoretical_spectrum
            .iter()
            .map(|mz| {
                let index = Self::calc_binned_position(*mz, bin_size, bin_offset);
                if preprocessed_experimental_spectrum[index] > 0.0 {
                    1
                } else {
                    0
                }
            })
            .sum()
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

        let ions_matched = Self::matched_ions(
            &theoretical_spectrum,
            &self.preprocessed_experimental_spectrum,
            self.config.bin_size,
            self.config.bin_offset,
        );

        // Score
        let score = Self::xcorr_spectra(
            &theoretical_spectrum,
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
    use std::{env, io::Write};

    use crate::utils::tests::{get_eng_experimental_spectrum, get_eng_fast_xcorr_spectrum};

    use itertools::multiunzip;
    use ndarray_stats::DeviationExt;
    use polars::prelude::*;
    use rayon::prelude::*;

    use crate::{
        configuration::{Configuration, FinalizedConfiguration},
        fast_xcorr::FastXcorr,
        utils::tests::{get_spectrum, read_test_data},
    };

    use super::*;

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

        let xcorr = FastXcorr::new(
            &config,
            (&experimental_spectrum.0, &experimental_spectrum.1),
            1,
        )
        .unwrap();

        let mut rouneded_xcorr_sped = xcorr
            .preprocessed_experimental_spectrum
            .iter()
            .map(|x| (x * 100.0).round() / 100.0)
            .collect::<Array1<f64>>();

        // Just select the values that are in the expected spectrum from Eng
        rouneded_xcorr_sped =
            rouneded_xcorr_sped.select(Axis(0), expected_xcorr_spec.0.as_slice().unwrap());

        let rounded_expected_xcorr_spec = expected_xcorr_spec
            .1
            .iter()
            .map(|x| (x * 100.0).round() / 100.0)
            .collect::<Array1<f64>>();

        if env::var("VERBOSE").is_ok() {
            let mut binned_experimental_spectrum = FastXcorr::experimental_spectrum_binning(
                &experimental_spectrum.0,
                &experimental_spectrum.1,
                config.bin_size,
                config.bin_offset,
                config.use_flanking_peaks,
            )
            .unwrap();
            binned_experimental_spectrum = binned_experimental_spectrum
                .select(Axis(0), expected_xcorr_spec.0.as_slice().unwrap());

            let output_file =
                std::fs::File::create("DIGSETK_fast_xcorr_preprocessing.tsv").unwrap();
            let mut output_writer = std::io::BufWriter::new(output_file);

            let _ = output_writer
                .write("bin\texpected\tcalc_binned\tcalc_prep\n".as_bytes())
                .unwrap();
            for (idx, (expected, (calc_binned, calc_prep))) in rounded_expected_xcorr_spec
                .iter()
                .zip(
                    binned_experimental_spectrum
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

        let config: FinalizedConfiguration = Configuration {
            bin_size: 1.0005,
            bin_offset: 0.4,
            ..Configuration::default()
        }
        .into();

        let xcorr = FastXcorr::new(
            &config,
            (&experimental_spectrum.0, &experimental_spectrum.1),
            1,
        )
        .unwrap();

        if env::var("VERBOSE").is_ok() {
            let peptide = CompoundPeptidoformIon::pro_forma("DIGSETK", None).unwrap();

            let mut binned_theoretical_spectrum =
                Array1::zeros(xcorr.preprocessed_experimental_spectrum.len());
            for mz in &xcorr.create_threoretical_spectrum(&peptide).unwrap() {
                let index =
                    FastXcorr::calc_binned_position(*mz, config.bin_size, config.bin_offset);
                binned_theoretical_spectrum[index] = 50.0;
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

        let scoring = xcorr.xcorr_peptide("DIGSETK").unwrap();
        println!("{scoring}");
        assert_eq!((scoring.score * 100.0).round() / 100.0, 2.92);
    }

    // Test xcorr implementations agains high-res MS data
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
                let charge = comet_df["charge"].i64().unwrap().get(idx).unwrap() as usize;

                let (mz_array, intensity_array) = get_spectrum(scan.to_string().as_str());

                let config: FinalizedConfiguration = Configuration {
                    use_flanking_peaks: true,
                    max_fragment_charge: 5,
                    ..Configuration::default()
                }
                .into();

                // fast xcorr implementation
                let fast_xcorr =
                    FastXcorr::new(&config, (&mz_array, &intensity_array), charge).unwrap();

                let scoring = fast_xcorr.xcorr_peptide(proforma_peptide).unwrap();

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

        if env::var("VERBOSE").is_ok() {
            let max_comet_xcorr = xcorrrs_df["comet_xcorr"].f64().unwrap().max().unwrap();

            let max_calculated_xcorr = xcorrrs_df["rs_xcorr"].f64().unwrap().max().unwrap();

            let mut plot = plotly::Plot::new();
            let diagonal_trace =
                plotly::Scatter::new(vec![0.0, max_comet_xcorr], vec![0.0, max_calculated_xcorr])
                    .mode(plotly::common::Mode::Lines)
                    .marker(plotly::common::Marker::default().color("red"))
                    .hover_info(plotly::common::HoverInfo::None)
                    .show_legend(false);

            let correlation_trace = plotly::Scatter::new(
                xcorrrs_df["comet_xcorr"].f64().unwrap().to_vec(),
                xcorrrs_df["rs_xcorr"].f64().unwrap().to_vec(),
            )
            .mode(plotly::common::Mode::Markers)
            .marker(plotly::common::Marker::default().color("blue"))
            .show_legend(false);

            plot.add_trace(diagonal_trace);
            plot.add_trace(correlation_trace);

            plot.set_layout(
                plotly::Layout::new()
                    .title("Comet xcorr vs rs_xcorr")
                    .x_axis(
                        plotly::layout::Axis::new()
                            .title("Comet xcorr")
                            .constrain(plotly::layout::AxisConstrain::Domain),
                    )
                    .y_axis(
                        plotly::layout::Axis::new()
                            .title("rs_xcorr")
                            .scale_anchor("x"),
                    ),
            );
            plot.write_html("99-fast_xcorrrs_vs_comet_xcorr.html");
        }

        // Normalize comet xcorrs and calculates xcorrs

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

        let scaled_comet_xcorr = xcorrrs_df
            .column("comet_xcorr")
            .unwrap()
            .f64()
            .unwrap()
            .to_ndarray()
            .unwrap()
            .mapv(|x| x / max_score);

        let scaled_fast_xcorrrs = xcorrrs_df
            .column("rs_xcorr")
            .unwrap()
            .f64()
            .unwrap()
            .to_ndarray()
            .unwrap()
            .mapv(|x| x / max_score);

        let rmse_fast_xcorrs = scaled_comet_xcorr
            .mean_sq_err(&scaled_fast_xcorrrs)
            .unwrap();

        println!("RMSE comet xcorr vs fast xcorrrs: {rmse_fast_xcorrs}");

        assert!(
            rmse_fast_xcorrs < 0.0002,
            "fast xcorr RMSE {rmse_fast_xcorrs} >= 0.0002"
        );
    }
}
