use ndarray::{s, Array1, Axis, Slice};
use ndrustfft::{ndfft, ndifft, Complex, FftHandler};
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

// pub trait IsXcorr {
//     fn config(&self) -> &Configuration;

//     fn max_experimental_mz(&self) -> f64;

//     fn prepared_experimental_spectrum(
//         &self,
//     ) -> &Array1<f64>;

//     fn xcorr_binned_spectrum(
//         binned_experimental_spectrum: &Array1<f64>,
//         binned_theoretical_spectrum: &Array1<f64>,
//     ) -> Result<f64, Error>;

//     /// Calculates the xcorr between a peptide and the experimental spectrum
//     ///
//     /// # Arguments
//     /// * `peptide` - The peptide sequence to score.
//     /// * `charge` - Charge states to use for scoring. If not provided, the configuration's charge states are used.
//     ///
//     fn xcorr_peptide(&self, peptide: &str, charge: usize) -> Result<ScoringResult, Error> {
//         let peptide = CompoundPeptidoformIon::pro_forma(peptide, None)
//             .map_err(Error::InvalidPeptideSequence)?;

//         let (min_theoretical_mass, max_theoretical_mass) =
//             match peptide.formulas().mass_bounds().into_option() {
//                 Some((min, max)) => (min.monoisotopic_mass().value, max.monoisotopic_mass().value),
//                 None => (-1.0, -1.0),
//             };

//         let mut fragment_charge = (charge - 1).max(1);
//         if fragment_charge > self.config().max_fragment_charge {
//             fragment_charge = self.config().max_fragment_charge;
//         }

//         let fragments = peptide.generate_theoretical_fragments(
//             Charge::new::<e>(fragment_charge),
//             &self.config().fragmentation_model,
//         );

//         let theoretical_spectrum = create_threoretical_spectrum(
//             &fragments,
//             fragment_charge,
//             self.max_experimental_mz(),
//             self.config().fragment_mass_range.clone(),
//         )?;
//         let ions_total = theoretical_spectrum.len();

//         let binned_thereoretical_spectrum = theoretical_spectrum_binning(
//             &theoretical_spectrum,
//             self.config().bin_size,
//             0, // Theoretical spectra are not shifted
//             Some(self.max_experimental_mz()),
//         )?;
//         drop(theoretical_spectrum);

//         // Score
//         let score = Self::xcorr_binned_spectrum(
//             self.prepared_experimental_spectrum(),
//             &binned_thereoretical_spectrum,
//         )?;

//         Ok(ScoringResult {
//             score,
//             min_theoretical_mass,
//             max_theoretical_mass,
//             ions_total,
//             ions_matched: 0,
//         })
//     }
// }

pub struct Xcorr<'a> {
    config: &'a Configuration,
    max_experimental_mz: f64,
    binned_experimental_spectrum: Array1<f64>,
}

impl Xcorr<'_> {
    /// Creates a new Xcorr instance.
    ///
    /// Arguments:
    /// * `config` - The configuration to use for scoring.
    /// * `experimental_spectrum` - The experimental spectrum to score against.
    ///
    pub fn new<'a>(
        config: &'a Configuration,
        experimental_spectrum: (&'a Array1<f64>, &'a Array1<f64>),
    ) -> Result<Xcorr<'a>, Error> {
        if experimental_spectrum.0.is_empty() {
            return Err(Error::EmptyExperimentalSpectrum);
        }

        if experimental_spectrum.0.len() != experimental_spectrum.1.len() {
            return Err(Error::ExperimentalSpectrumShape(
                experimental_spectrum.0.len(),
                experimental_spectrum.1.len(),
            ));
        }

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

        let max_experimental_mz = filtered_experimental_spectrum
            .0
            .iter()
            .fold(f64::NEG_INFINITY, |a, &b| a.max(b));

        let shift = (MZ_SHIFT as f64 / config.bin_size) as usize;
        let binned_experimental_spectrum = experimental_spectrum_binning(
            &filtered_experimental_spectrum.0,
            &filtered_experimental_spectrum.1,
            config.bin_size,
            shift,
        )?;

        Ok(Xcorr {
            config,
            max_experimental_mz,
            binned_experimental_spectrum,
        })
    }

    /// FTT based cross-correlation for binned spectra.
    /// Mostly copied from `fftconvolve`-crate.
    ///
    /// # Arguments
    /// * `binned_experimental_spectrum` - Binned experimental spectrum +/- 75mz shift
    /// * `binned_theoretical_spectrum` - Binned theoretical spectrum
    ///
    pub fn binned_spectrum_cross_correlation(
        binned_experimental_spectrum: &Array1<f64>,
        binned_theoretical_spectrum: &Array1<f64>,
    ) -> Array1<f64> {
        let padding = binned_experimental_spectrum.len() - binned_theoretical_spectrum.len();

        let mut experimental_fft =
            Array1::<Complex<f64>>::zeros(binned_experimental_spectrum.len());
        let experimental_complex = binned_experimental_spectrum.mapv(|v| Complex::new(v, 0.0));
        let handler: FftHandler<f64> = FftHandler::new(binned_experimental_spectrum.len());
        ndfft(&experimental_complex, &mut experimental_fft, &handler, 0);
        drop(experimental_complex);

        // Create complex array for the theoretical spectrum with padding
        let mut theoretical_complex =
            Array1::<Complex<f64>>::zeros(binned_theoretical_spectrum.len() + padding);
        theoretical_complex
            .slice_mut(s![0..binned_theoretical_spectrum.len()])
            .iter_mut()
            .zip(binned_theoretical_spectrum.iter())
            .for_each(|(c, &v)| *c = Complex::new(v, 0.0));

        let mut theoretical_fft = Array1::<Complex<f64>>::zeros(theoretical_complex.len());
        let handler: FftHandler<f64> = FftHandler::new(theoretical_complex.len());
        ndfft(&theoretical_complex, &mut theoretical_fft, &handler, 0);
        drop(theoretical_complex);

        theoretical_fft.mapv_inplace(|elem| elem.conj());
        let corr = experimental_fft * theoretical_fft;

        let mut corr_ifft = Array1::<Complex<f64>>::zeros(corr.len());
        let handler: FftHandler<f64> = FftHandler::new(corr_ifft.len());
        ndifft(&corr, &mut corr_ifft, &handler, 0);
        drop(corr);
        let mut corr = corr_ifft.into_iter().map(|c| c.re).collect::<Array1<f64>>();

        // Mode "valid"
        let new_shape = binned_experimental_spectrum.len() - binned_theoretical_spectrum.len() + 1;
        corr.slice_axis_inplace(Axis(0), Slice::new(0, Some(new_shape as isize), 1));
        corr
    }

    /// Calculates the xcorr between an already binned theoretical and binned experimental spectrum.
    ///
    /// # Arguments
    /// * `binned_experimental_spectrum` - The binned experimental m/z values + 75 m/z shift.
    /// * `binned_theoretical_spectrum` - The binned theoretical m/z value
    ///
    pub fn xcorr_binned_spectrum(
        binned_experimental_spectrum: &Array1<f64>,
        binned_theoretical_spectrum: &Array1<f64>,
    ) -> Result<f64, Error> {
        let correlation = Self::binned_spectrum_cross_correlation(
            binned_experimental_spectrum,
            binned_theoretical_spectrum,
        );

        if correlation.is_empty() {
            return Ok(-99999.0);
        }

        let middle = correlation.len() / 2;
        let zeroshift_corr = correlation[middle];
        let correlation = ndarray::concatenate(
            Axis(0),
            &[
                correlation.slice(s![0..middle]),
                correlation.slice(s![middle + 1..]),
            ],
        )
        .unwrap();

        // mean_corr = np.mean(corr) #Background similarity
        let mean_corr = match correlation.mean() {
            Some(mean) => mean,
            None => return Ok(-100000.0),
        };

        // xcorr_score = np.round((zeroshift_corr - mean_corr) / 10000, 4)
        Ok((zeroshift_corr - mean_corr) / 10000.0)
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
            0, // Theoretical spectra are not shifted
            Some(self.max_experimental_mz),
        )?;
        drop(theoretical_spectrum);

        // Score
        let score = Self::xcorr_binned_spectrum(
            &self.binned_experimental_spectrum,
            &binned_thereoretical_spectrum,
        )?;

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

    use ndarray_stats::DeviationExt;
    use polars::prelude::*;
    use rayon::prelude::*;

    use super::*;

    /// Test `binned_spectrum_cross_correlation` function with a veeeeeeeery simple example validated against numpy.
    ///
    #[test]
    fn test_binned_spectrum_cross_correlation() {
        let x: Array1<f64> = Array1::from(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let y: Array1<f64> = Array1::from(vec![1.0, 2.0, 3.0]);
        let result = Xcorr::binned_spectrum_cross_correlation(&x, &y);
        assert_eq!(result, Array1::from(vec![14.0, 20.0, 26.0, 32.0]));
    }

    // Test the xcorr calculation against Comet's xcorr values.
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
                let xcorr = Xcorr::new(&config, (&mz_array, &intensity_array)).unwrap();

                let scoring = match xcorr.xcorr_peptide(proforma_peptide, charge) {
                    Ok(scoring) => scoring,
                    Err(err) => {
                        println!(
                            "Empty theoretical spectrum for scan {} and peptide {}",
                            scan, err
                        );
                        panic!(
                            "Empty theoretical spectrum for scan {} and peptide {}",
                            scan, err
                        );
                    }
                };

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

        xcorrrs_df
            .sort_in_place(
                ["comet_xcorr"],
                SortMultipleOptions::default().with_order_descending(true),
            )
            .unwrap();

        CsvWriter::new(std::fs::File::create("xcorr_results.tsv").unwrap())
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
        plot.write_html("99-xcorr-correlation.html");

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
}
