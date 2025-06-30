use ndarray::{s, Array1, Axis};
use rustyms::{
    system::{e, usize::Charge},
    CompoundPeptidoformIon,
};

use crate::{
    binning::{experimental_spectrum_binning, theoretical_spectrum_binning},
    configuration::Configuration,
    error::Error,
    scoring_result::ScoringResult,
    utils::create_threoretical_spectrum,
};

/// +/- m/z shift for the y' calculation.
const MZ_SHIFT: u8 = 75;

pub struct FastXcorr<'a> {
    config: &'a Configuration,
    max_experimental_mz: f64,
    shift: usize,
    y_prime: Array1<f64>,
}

impl FastXcorr<'_> {
    /// Creates a new FastXcorr instance.
    ///
    /// Arguments:
    /// * `config` - The configuration to use for scoring.
    /// * `experimental_spectrum` - The experimental spectrum to score against.
    ///
    pub fn new<'a>(
        config: &'a Configuration,
        experimental_spectrum: (&'a Array1<f64>, &'a Array1<f64>),
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

        let shift = (MZ_SHIFT as f64 / config.bin_size) as usize;
        let binned_experimental_spectrum = experimental_spectrum_binning(
            &filtered_experimental_spectrum.0,
            &filtered_experimental_spectrum.1,
            config.bin_size,
            shift,
        )?;

        let y_prime = Self::calc_y_prime(&binned_experimental_spectrum, shift);

        Ok(FastXcorr {
            config,
            max_experimental_mz,
            shift,
            y_prime,
        })
    }

    /// Calculates the specturm +/- 75 m/z shift from equation 6 in https://pubs.acs.org/doi/10.1021/pr800420s
    ///
    /// Arguments:
    /// * `binned_experimental_spectrum` - The binned experimental m/z values +/- m/z shift.
    /// * `shift` - Number of bins to shift.
    ///
    pub fn calc_y_prime_shift(
        binned_experimental_spectrum: &Array1<f64>,
        shift: usize,
    ) -> Array1<f64> {
        // Extend by `shift` bins on both sides
        let mut y_prime_shift = Array1::zeros(binned_experimental_spectrum.len());

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
    /// * `shift` - Number of bins to shift.
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
        // NOTE: / 10000.0 is not mentioned in the paper, but seems to be necessary to match the results
        // Also the correlations-version is using this scaling
        binned_theoretical_spectrum.dot(y_prime) / 10000.0
    }

    /// Calculates the xcorr between a peptide and the experimental spectrum
    ///
    /// # Arguments
    /// * `peptide` - The peptide sequence to score.
    /// * `charge` - Charge states to use for scoring. If not provided, the configuration's charge states are used.
    ///
    pub fn xcorr_peptide(&self, peptide: &str, charge: usize) -> Result<ScoringResult, Error> {
        let peptide = CompoundPeptidoformIon::pro_forma(peptide, None)
            .map_err(Error::InvalidPeptideSequence)?;

        let (min_theoretical_mass, max_theoretical_mass) =
            match peptide.formulas().mass_bounds().into_option() {
                Some((min, max)) => (min.monoisotopic_mass().value, max.monoisotopic_mass().value),
                None => (-1.0, -1.0),
            };

        let mut fragment_charge = (charge - 1).max(1);
        if fragment_charge > self.config.max_fragment_charge {
            fragment_charge = self.config.max_fragment_charge;
        }

        let fragments = peptide.generate_theoretical_fragments(
            Charge::new::<e>(fragment_charge),
            &self.config.fragmentation_model,
        );

        let theoretical_spectrum =
            create_threoretical_spectrum(&fragments, fragment_charge, self.max_experimental_mz)?;
        let ions_total = theoretical_spectrum.len();

        let binned_thereoretical_spectrum = theoretical_spectrum_binning(
            &theoretical_spectrum,
            self.config.bin_size,
            self.shift,
            Some(self.max_experimental_mz),
        )?;

        drop(theoretical_spectrum);

        // Score
        let score = Self::xcorr_binned_spectrum(&binned_thereoretical_spectrum, &self.y_prime);

        Ok(ScoringResult {
            score,
            min_theoretical_mass,
            max_theoretical_mass,
            ions_total,
            ions_matched: 0,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::utils::tests::{get_spectrum, read_test_data};
    use std::{ops::Neg, path::PathBuf};

    use dihardts_omicstools::proteomics::io::mzml::reader::Reader as MzmlReader;
    use ndarray_stats::DeviationExt;
    use polars::prelude::*;
    use rayon::prelude::*;

    use super::*;

    #[test]
    fn test_y_prime() {
        let shift: usize = 5;
        let exp = ndarray::concatenate(
            Axis(0),
            &[
                Array1::from(vec![0.0; shift]).view(),
                Array1::from(vec![1.0; 16]).view(),
                Array1::from(vec![0.0; shift]).view(),
            ],
        )
        .unwrap();

        let y_prime_shift = FastXcorr::calc_y_prime_shift(&exp, shift);

        let expected_y_shift = Array1::from(vec![
            1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
            9.0, 8.0, 7.0, 6.0, 5.0, 5.0, 4.0, 3.0, 2.0, 1.0,
        ]);

        assert_eq!(y_prime_shift, expected_y_shift);

        let y_prime = FastXcorr::calc_y_prime(&exp, shift);

        let expected_y_prime = Array1::from(vec![
            -0.006666666666666667,
            -0.013333333333333334,
            -0.02,
            -0.02666666666666667,
            -0.03333333333333333,
            0.9666666666666667,
            0.96,
            0.9533333333333334,
            0.9466666666666667,
            0.94,
            0.9333333333333333,
            0.9333333333333333,
            0.9333333333333333,
            0.9333333333333333,
            0.9333333333333333,
            0.9333333333333333,
            0.94,
            0.9466666666666667,
            0.9533333333333334,
            0.96,
            0.9666666666666667,
            -0.03333333333333333,
            -0.02666666666666667,
            -0.02,
            -0.013333333333333334,
            -0.006666666666666667,
        ]);

        assert_eq!(y_prime, expected_y_prime);
    }

    /// Test the xcorr calculation against Comet's xcorr values.
    ///
    #[test]
    fn test_xcorr() {
        let comet_df = read_test_data();

        let number_of_psms = match std::env::var("TEST_NUMBER_OF_PSMS") {
            Ok(env_value) => env_value.parse::<usize>().unwrap(),
            Err(_) => comet_df.height(),
        };

        #[allow(clippy::type_complexity)]
        let (scan_col, (peptide_col, (comet_xcorr_col, xcorrrs_col))): (
            Vec<i64>,
            (Vec<String>, (Vec<f64>, Vec<f64>)),
        ) = (0..number_of_psms)
            .into_par_iter()
            .map(|idx| {
                let scan = comet_df["scan"].i64().unwrap().get(idx).unwrap();
                let comet_xcorr = comet_df["xcorr"].f64().unwrap().get(idx).unwrap();
                let proforma_peptide = comet_df["proforma_peptide"]
                    .str()
                    .unwrap()
                    .get(idx)
                    .unwrap();
                let charge = comet_df["charge"].i64().unwrap().get(idx).unwrap() as usize;

                let (mz_array, intensity_array) = get_spectrum(scan.to_string().as_str());

                let config = Configuration::default();
                let xcorr = FastXcorr::new(&config, (&mz_array, &intensity_array)).unwrap();

                let scoring = xcorr.xcorr_peptide(proforma_peptide, charge).unwrap();

                (
                    scan,
                    (proforma_peptide.to_string(), (comet_xcorr, scoring.score)),
                )
            })
            .unzip();

        let mut xcorrrs_df = DataFrame::new(vec![
            Column::new("scan".into(), scan_col),
            Column::new("modified_peptide".into(), peptide_col),
            Column::new("comet_xcorr".into(), comet_xcorr_col),
            Column::new("xcorrrs".into(), xcorrrs_col),
        ])
        .unwrap();

        CsvWriter::new(std::fs::File::create("fast_xcorr_results.tsv").unwrap())
            .with_separator(b'\t')
            .finish(&mut xcorrrs_df)
            .unwrap();

        // Plot the comet xcorrs vs xcorrrs
        let max_comet_xcorr = xcorrrs_df["comet_xcorr"].f64().unwrap().max().unwrap();
        let max_calculated_xcorr = xcorrrs_df["xcorrrs"].f64().unwrap().max().unwrap();
        let max_score = max_comet_xcorr.max(max_calculated_xcorr);

        let mut plot = plotly::Plot::new();
        let diagonal_trace = plotly::Scatter::new(vec![0.0, max_score], vec![0.0, max_score])
            .mode(plotly::common::Mode::Lines)
            .marker(plotly::common::Marker::default().color("red"))
            .hover_info(plotly::common::HoverInfo::None)
            .show_legend(false);

        let correlation_trace = plotly::Scatter::new(
            xcorrrs_df["comet_xcorr"].f64().unwrap().to_vec(),
            xcorrrs_df["xcorrrs"].f64().unwrap().to_vec(),
        )
        .mode(plotly::common::Mode::Markers)
        .marker(plotly::common::Marker::default().color("blue"))
        .show_legend(false);

        plot.add_trace(diagonal_trace);
        plot.add_trace(correlation_trace);

        plot.set_layout(
            plotly::Layout::new()
                .title("Comet xcorr vs xcorrrs")
                .x_axis(
                    plotly::layout::Axis::new()
                        .title("Comet xcorr")
                        .constrain(plotly::layout::AxisConstrain::Domain),
                )
                .y_axis(
                    plotly::layout::Axis::new()
                        .title("xcorrrs")
                        .scale_anchor("x"),
                ),
        );
        plot.write_html("99-fast-xcorr-correlation.html");

        // Normalize comet xcorrs and calculates xcorrs

        let comet_xcorr_max = xcorrrs_df
            .column("comet_xcorr")
            .unwrap()
            .f64()
            .unwrap()
            .max()
            .unwrap();

        let xcorrrs_max = xcorrrs_df
            .column("xcorrrs")
            .unwrap()
            .f64()
            .unwrap()
            .max()
            .unwrap();

        let max_score = comet_xcorr_max.max(xcorrrs_max);

        let scaled_comet_xcorr = xcorrrs_df
            .column("comet_xcorr")
            .unwrap()
            .f64()
            .unwrap()
            .to_ndarray()
            .unwrap()
            .mapv(|x| x / max_score);

        let scaled_xcorrrs = xcorrrs_df
            .column("xcorrrs")
            .unwrap()
            .f64()
            .unwrap()
            .to_ndarray()
            .unwrap()
            .mapv(|x| x / max_score);

        let rmse = scaled_comet_xcorr
            .root_mean_sq_err(&scaled_xcorrrs)
            .unwrap();

        // RMSE under 0.02 should be good enough
        assert!(rmse < 0.02, "RMSE > 0.02: {}", rmse);
    }

    /// Creates plots of the the spectra for visual inspection.
    /// Skipped by default
    #[test]
    #[ignore = "Visual inspection test, skipped by default"]
    fn test_prints() {
        let prints_folder = PathBuf::from("./prints");
        if prints_folder.exists() {
            std::fs::remove_dir_all(&prints_folder).unwrap();
        }
        std::fs::create_dir(&prints_folder).unwrap();

        let comet_df = read_test_data();

        let mut mzml_byte_reader = std::io::BufReader::new(
            std::fs::File::open("test_files/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.mzML")
                .unwrap(),
        );
        let mut mzml = MzmlReader::read_indexed(&mut mzml_byte_reader, None, true, false).unwrap();

        let config = Configuration::default();

        let number_of_psms = match std::env::var("TEST_NUMBER_OF_PSMS") {
            Ok(env_value) => env_value.parse::<usize>().unwrap(),
            Err(_) => comet_df.height(),
        };

        for i in 0..number_of_psms {
            let scan = comet_df["scan"].i64().unwrap().get(i).unwrap();
            let _xcorr = comet_df["xcorr"].f64().unwrap().get(i).unwrap();
            let proforma_peptide = comet_df["proforma_peptide"].str().unwrap().get(i).unwrap();
            let charge = comet_df["charge"].i64().unwrap().get(i).unwrap() as usize;
            let exp_neutral_mass = comet_df["exp_neutral_mass"].f64().unwrap().get(i).unwrap();

            let plot_base_path = prints_folder.join(format!("{}_scan_{}", i, scan));
            std::fs::create_dir_all(&plot_base_path).unwrap();

            let peptide = CompoundPeptidoformIon::pro_forma(proforma_peptide, None).unwrap();

            let mut fragment_charge = (charge - 1).max(1);
            if fragment_charge > config.max_fragment_charge {
                fragment_charge = config.max_fragment_charge;
            }

            let fragments = peptide.generate_theoretical_fragments(
                Charge::new::<e>(fragment_charge),
                &config.fragmentation_model,
            );

            // figure 1a
            let binary_data_array_list = mzml
                .get_spectrum(&format!(
                    "controllerType=0 controllerNumber=1 scan={}",
                    scan
                ))
                .unwrap()
                .binary_data_array_list;

            let mz_array = Array1::from(
                binary_data_array_list
                    .get_mz_array()
                    .unwrap()
                    .deflate_data()
                    .unwrap(),
            );
            let intensity_array = Array1::from(
                binary_data_array_list
                    .get_intensity_array()
                    .unwrap()
                    .deflate_data()
                    .unwrap(),
            );

            let max_intensity = intensity_array
                .iter()
                .fold(f64::NEG_INFINITY, |a, &b| a.max(b));

            let max_mz = mz_array.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

            let theoretical_spectrum =
                create_threoretical_spectrum(&fragments, fragment_charge, max_mz).unwrap();

            let experimental_trace = plotly::Bar::new(mz_array.to_vec(), intensity_array.to_vec())
                .name("orignal spectrum")
                .width(2.0);

            let theoretical_trace = plotly::Bar::new(
                theoretical_spectrum.to_vec(),
                vec![max_intensity.neg(); theoretical_spectrum.len()],
            )
            .name("theoretical spectrum")
            .width(2.0);

            let mut plot: plotly::Plot = plotly::Plot::new();
            plot.add_trace(experimental_trace);
            plot.add_trace(theoretical_trace);
            plot.set_layout(
                plotly::Layout::new()
                    .title(format!(
                        "{} / {} Da / {} charge - original spectrum",
                        proforma_peptide, exp_neutral_mass, charge,
                    ))
                    .box_gap(0.4),
            );
            plot.write_html(plot_base_path.join("1a-original-spectrum.html"));

            // figure 1b
            let experimental_trace =
                plotly::Bar::new(mz_array.to_vec(), intensity_array.sqrt().to_vec())
                    .name("orignal spectrum")
                    .width(2.0);

            let theoretical_trace = plotly::Bar::new(
                theoretical_spectrum.to_vec(),
                vec![max_intensity.sqrt().neg(); theoretical_spectrum.len()],
            )
            .name("theoretical spectrum")
            .width(2.0);

            let mut plot: plotly::Plot = plotly::Plot::new();
            plot.add_trace(experimental_trace);
            plot.add_trace(theoretical_trace);
            plot.set_layout(
                plotly::Layout::new()
                    .title(format!(
                        "{} / {} Da / {} charge - square root intensities",
                        proforma_peptide, exp_neutral_mass, charge,
                    ))
                    .box_gap(0.4),
            );
            plot.write_html(plot_base_path.join("1b-sqrt-intensities.html"));

            // figure 1c
            let binned_experimental_spectrum =
                experimental_spectrum_binning(&mz_array, &intensity_array, 0.02, 0).unwrap();

            let binned_theoretical_spectrum = theoretical_spectrum_binning(
                &theoretical_spectrum,
                0.02,
                0,
                Some(exp_neutral_mass),
            )
            .unwrap()
            .neg();

            let experimental_trace = plotly::Bar::new(
                (0..binned_experimental_spectrum.len() as u32)
                    .map(|i| i as f64 * 0.02)
                    .collect::<Vec<f64>>(),
                binned_experimental_spectrum.to_vec(),
            );

            let theoretical_trace = plotly::Bar::new(
                (0..binned_theoretical_spectrum.len() as u32)
                    .map(|i| i as f64 * 0.02)
                    .collect::<Vec<f64>>(),
                binned_theoretical_spectrum.to_vec(),
            );

            let mut plot: plotly::Plot = plotly::Plot::new();
            plot.add_trace(experimental_trace);
            plot.add_trace(theoretical_trace);
            plot.set_layout(
                plotly::Layout::new()
                    .title(format!(
                        "{} / {} Da / {} charge - intensities normalized spectrum",
                        proforma_peptide, exp_neutral_mass, charge,
                    ))
                    .box_gap(0.4),
            );
            plot.write_html(plot_base_path.join("1c-intensities-normalized-spectrum.html"));

            // figure 1d
            let xcorr = FastXcorr::new(&config, (&mz_array, &intensity_array)).unwrap();

            let fast_xcorr_spectrum = &xcorr.y_prime;

            let shifte_accomodated_theoretical_spectrum = ndarray::concatenate(
                Axis(0),
                &[
                    Array1::zeros(xcorr.shift).view(),
                    binned_theoretical_spectrum.view(),
                    Array1::zeros(xcorr.shift).view(),
                ],
            )
            .unwrap();

            let experimental_trace = plotly::Bar::new(
                (0..fast_xcorr_spectrum.len() as u32)
                    .map(|i| i as f64 * 0.02)
                    .collect::<Vec<f64>>(),
                fast_xcorr_spectrum.to_vec(),
            );

            let theoretical_trace = plotly::Bar::new(
                (0..shifte_accomodated_theoretical_spectrum.len() as u32)
                    .map(|i| i as f64 * 0.02)
                    .collect::<Vec<f64>>(),
                shifte_accomodated_theoretical_spectrum.to_vec(),
            );

            let mut plot: plotly::Plot = plotly::Plot::new();
            plot.add_trace(experimental_trace);
            plot.add_trace(theoretical_trace);
            plot.set_layout(
                plotly::Layout::new()
                    .title(format!(
                        "{} / {} Da / {} charge - fast xcorr spectrum",
                        proforma_peptide, exp_neutral_mass, charge,
                    ))
                    .box_gap(0.4),
            );
            plot.write_html(plot_base_path.join("1d-fast-xcorr-spectrum.html"));
        }
    }
}
