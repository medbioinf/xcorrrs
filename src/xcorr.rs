use ndarray::{s, Array1, Axis, Slice};
use ndrustfft::{ndfft, ndifft, Complex, FftHandler};
use rustyms::CompoundPeptidoformIon;

use crate::{
    binning::{experimental_spectrum_binning, theoretical_spectrum_binning},
    configuration::Configuration,
    error::Error,
    scoring_result::ScoringResult,
    utils::{calculate_number_of_bins_to_shift, create_threoretical_spectrum},
};

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
    shift: usize,
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
        charge: usize,
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

        Ok(Xcorr {
            config,
            max_experimental_mz,
            shift,
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

        let theoretical_spectrum = create_threoretical_spectrum(
            &peptide,
            &self.config.fragmentation_model,
            fragment_charge,
            self.max_experimental_mz,
        )?;

        let ions_total = theoretical_spectrum.len();

        let binned_thereoretical_spectrum = theoretical_spectrum_binning(
            &theoretical_spectrum,
            self.config.bin_size,
            self.config.bin_offset,
            charge,
            0, // Theoretical spectra are not shifted
            Some(self.max_experimental_mz),
        )?;
        drop(theoretical_spectrum);

        let ions_matched = self
            .binned_experimental_spectrum
            .slice(s![
                self.shift..self.binned_experimental_spectrum.len() - self.shift
            ])
            .iter()
            .zip(binned_thereoretical_spectrum.iter())
            .fold(0_usize, |count, (experimental_peak, theoretical_peak)| {
                if *experimental_peak <= 0.0 || *theoretical_peak <= 0.0 {
                    count
                } else {
                    count + 1
                }
            });

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
            ions_matched,
        })
    }
}
