use std::fmt::Display;

use context_error::{BasicKind, BoxedError};
use mzannotate::{
    annotation::model::{FragmentationModel, PrimaryIonSeries},
    fragment::{Fragment, FragmentType},
    prelude::PeptidoformFragmentation,
};
use mzcore::{
    chemistry::{Element::H as Hydrogen, MassMode, NeutralLoss, OutputMolecularFormula},
    molecular_formula,
    sequence::PeptidoformIon,
    system::{e, isize::Charge},
};
use mzdata::prelude::{ByteArrayView, IonProperties, PrecursorSelection, SpectrumLike};
use mzdata::spectrum::{ArrayType, SignalContinuity};
use mzpeaks::{CentroidLike, DeconvolutedCentroidLike};
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Resolves the fragment/precursor charge to use for a spectrum: the spectrum's own precursor
/// charge takes precedence, falling back to `charge_fallback` only when the spectrum doesn't
/// report one (e.g. it wasn't determined by the instrument/file).
pub fn resolve_precursor_charge<C, D, S>(spectrum: &S, charge_fallback: Option<usize>) -> Option<usize>
where
    C: CentroidLike,
    D: DeconvolutedCentroidLike,
    S: SpectrumLike<C, D>,
{
    spectrum
        .precursor()
        .and_then(|precursor| precursor.ion())
        .and_then(|ion| ion.charge())
        .and_then(|charge| usize::try_from(charge).ok())
        .or(charge_fallback)
}

/// Marker error returned by [`extract_peaks`] when the spectrum's signal continuity is
/// `Profile` - i.e. it hasn't been centroided/peak-picked yet. Both `xcorr` and `msamanda`
/// assume a discrete peak list; treating a continuous profile trace as one would silently
/// count every raw data point as its own "peak", corrupting the score. Each caller maps this
/// to its own `Error::UncentroidedSpectrum` variant.
#[derive(Debug)]
pub struct UncentroidedSpectrum;

/// Extracts `(m/z, intensity)` peaks from a spectrum, preserving full `f64` intensity precision
/// whenever possible.
///
/// Returns [`UncentroidedSpectrum`] if `spectrum.signal_continuity()` is `Profile`; `Unknown`
/// and `Centroid` are both accepted, since not every centroided source (e.g. this crate's own
/// [`crate::utils::tests::build_raw_spectrum`] test fixtures) bothers to tag that metadata.
///
/// `mzpeaks`'s centroided peak types (and the peak-agnostic [`mzdata`] peak iterator built on
/// top of them) store intensity as `f32`, which is fine once a spectrum has actually been
/// centroided (that's the precision it was picked at), but needlessly truncates a spectrum that
/// still carries its raw/profile data arrays, which are commonly `f64`. So raw arrays are
/// preferred here when present; only spectra that have already been reduced to centroided peaks
/// (no raw arrays left) fall back to the `f32` peak iterator.
pub fn extract_peaks<C, D, S>(spectrum: &S) -> Result<(Vec<f64>, Vec<f64>), UncentroidedSpectrum>
where
    C: CentroidLike,
    D: DeconvolutedCentroidLike,
    S: SpectrumLike<C, D>,
{
    if spectrum.signal_continuity() == SignalContinuity::Profile {
        return Err(UncentroidedSpectrum);
    }

    let from_raw_arrays = spectrum.raw_arrays().and_then(|arrays| {
        let mz = arrays.mzs().ok()?;
        let intensity = arrays.get(&ArrayType::IntensityArray)?.to_f64().ok()?;
        Some((mz.into_owned(), intensity.into_owned()))
    });

    Ok(from_raw_arrays.unwrap_or_else(|| {
        spectrum
            .peaks()
            .iter()
            .map(|peak| (peak.mz, peak.intensity as f64))
            .unzip()
    }))
}

#[derive(Error, Debug)]
pub enum Error {
    #[error("Invalid ProForma string: {0:?}")]
    InvalidProForma(Vec<BoxedError<'static, BasicKind>>),
    #[error("Invalid configuration: {0}")]
    InvalidConfiguration(String),
}

/// Ion types used in configuration to select which ion series to use
/// for the fragmentation model. Shared between the `xcorr` and `msamanda` scorers.
///
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum Ion {
    A,
    B,
    C,
    X,
    Y,
    Z,
}

impl Display for Ion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Ion::A => write!(f, "A"),
            Ion::B => write!(f, "B"),
            Ion::C => write!(f, "C"),
            Ion::X => write!(f, "X"),
            Ion::Y => write!(f, "Y"),
            Ion::Z => write!(f, "Z"),
        }
    }
}

impl TryFrom<String> for Ion {
    type Error = Error;

    fn try_from(value: String) -> Result<Self, Self::Error> {
        match value.to_ascii_uppercase().as_str() {
            "A" => Ok(Ion::A),
            "B" => Ok(Ion::B),
            "C" => Ok(Ion::C),
            "X" => Ok(Ion::X),
            "Y" => Ok(Ion::Y),
            "Z" => Ok(Ion::Z),
            _ => Err(Self::Error::InvalidConfiguration(format!(
                "Unknown ion type: {value}",
            ))),
        }
    }
}

/// Builds a `FragmentationModel` from the selected ion series and whether to include
/// the common neutral-loss variants (water/ammonia) on b/y ions.
///
pub fn build_fragmentation_model(ions: &[Ion], use_neutral_loss_ions: bool) -> FragmentationModel {
    let mut fragmentation_model = FragmentationModel::none().clone();

    if ions.contains(&Ion::A) {
        fragmentation_model.a = PrimaryIonSeries::default();
    };

    if ions.contains(&Ion::B) {
        fragmentation_model.b = PrimaryIonSeries::default();
    };

    if ions.contains(&Ion::C) {
        fragmentation_model.c = PrimaryIonSeries::default();
    };

    if ions.contains(&Ion::X) {
        fragmentation_model.x = PrimaryIonSeries::default();
    };

    if ions.contains(&Ion::Y) {
        fragmentation_model.y = PrimaryIonSeries::default();
    };

    if ions.contains(&Ion::Z) {
        fragmentation_model.z = PrimaryIonSeries::default();
    };

    if use_neutral_loss_ions {
        let mut b_ions = fragmentation_model.b.clone();
        b_ions = b_ions.neutral_losses(vec![
            NeutralLoss::Loss(1, molecular_formula!(H 2 O 1)),
            NeutralLoss::Loss(1, molecular_formula!(N 1 H 3)),
        ]);
        fragmentation_model.b = b_ions;

        let mut y_ions = fragmentation_model.y.clone();
        y_ions = y_ions.neutral_losses(vec![
            NeutralLoss::Loss(1, molecular_formula!(H 2 O 1)),
            NeutralLoss::Loss(1, molecular_formula!(N 1 H 3)),
        ]);
        fragmentation_model.b = y_ions
    }

    fragmentation_model
}

/// Da to m/z conversion
///
/// # Arguments
/// * `mass` - Mass in Dalton
/// * `charge` - Charge state
///
pub fn dalton_to_mass_to_charge(mass: f64, charge: usize) -> f64 {
    let charge = charge as f64;
    (mass + Hydrogen.mass(None).unwrap().value * charge) / charge
}

/// Converts mass to charge ration (Thompson) to Dalton
///
/// # Arguments
/// * `mz` - Mass to charge ratio (Thompson)
/// * `charge` - Charge
///
pub fn mass_to_charge_to_dalton(mz: f64, charge: usize) -> f64 {
    let charge = charge as f64;
    mz * charge - Hydrogen.mass(None).unwrap().value * charge
}

/// Creates theoretical fragments for the given peptide using the provided fragmentation model and maximum charge.
///
/// # Arguments
/// * `peptide` - The peptide for which to generate the theoretical spectrum.
/// * `fragmentation_model` - The fragmentation model to use for generating fragments.
/// * `max_charge` - The maximum charge state to consider for the fragments.
///
pub fn create_theoretical_fragments(
    peptidoform_ion: &PeptidoformIon,
    fragmentation_model: &FragmentationModel,
    max_charge: usize,
) -> Result<Vec<Fragment<OutputMolecularFormula>>, Error> {
    let mut fragments: Vec<Fragment<OutputMolecularFormula>> = peptidoform_ion
        .generate_theoretical_fragments(Charge::new::<e>(max_charge as isize), fragmentation_model)
        .into_iter()
        .filter_map(|f| {
            // `mzannotate` unconditionally appends one intact-precursor "fragment" per
            // peptide regardless of the fragmentation model; Comet's xcorr ion ladder
            // never includes such a peak, so it must be excluded here.
            if matches!(f.ion, FragmentType::Precursor) {
                None
            } else if f.mz(MassMode::Monoisotopic).is_some() {
                Some(f)
            } else {
                None
            }
        })
        .collect();

    fragments.sort_by(|a, b| {
        a.mz(MassMode::Monoisotopic)
            .partial_cmp(&b.mz(MassMode::Monoisotopic))
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    Ok(fragments)
}

#[cfg(test)]
pub mod tests {
    use crate::xcorr::configuration::{Configuration, FinalizedConfiguration};

    use mzcore::ontology::Ontologies;

    use super::*;

    use mzdata::spectrum::{
        ArrayType, BinaryArrayMap, BinaryDataArrayType, DataArray, Precursor, RawSpectrum,
        SelectedIon, SpectrumDescription,
    };

    /// Just a sanity check to make sure that max charge is used as is in mzannotate
    #[test]
    fn test_fragment_creation() {
        let max_charge = 6;
        let config: FinalizedConfiguration = Configuration::default().into();
        let peptide = PeptidoformIon::pro_forma("DIGSETK", &Ontologies::empty())
            .unwrap()
            .0;

        let mut charge_states = vec![false; max_charge + 1];
        charge_states[0] = true; // Charge 0 is not valid, but we use it for indexing

        peptide
            .generate_theoretical_fragments::<OutputMolecularFormula>(
                Charge::new::<e>(max_charge as isize),
                &config.fragmentation_model,
            )
            .iter()
            .for_each(|f| charge_states[f.charge.value as usize] = true);

        assert!(charge_states.iter().all(|x| *x));
    }

    /// Reads a spectrum from its postcard/zstd fixture file.
    ///
    /// # Arguments
    /// * `scan` - The scan number to read the spectrum for.
    ///
    pub fn get_spectrum(scan: &str) -> mzdata::spectrum::Spectrum {
        let compressed =
            std::fs::read(format!("test_files/spectra/scan_{scan}.postcard.zst")).unwrap();
        let bytes = zstd::decode_all(compressed.as_slice()).unwrap();
        postcard::from_bytes(&bytes).unwrap()
    }

    /// Builds a raw (unpicked) [`mzdata`] spectrum from parallel m/z + intensity arrays, plus the
    /// precursor charge and m/z, for feeding into [`crate::xcorr::fast_xcorr::FastXcorr::new`] /
    /// [`crate::msamanda::deconvolution::Spectrum::from_mzdata_spectrum`] in tests.
    ///
    /// Stores both arrays as `f64` `DataArray`s (not `mzpeaks::CentroidPeak`, which is `f32`
    /// intensity only) so tests retain the full double-precision intensity that
    /// `crate::utils::extract_peaks` prefers when raw arrays are present.
    pub fn build_raw_spectrum(
        mz: &[f64],
        intensity: &[f64],
        precursor_charge: i32,
        precursor_mz: f64,
    ) -> RawSpectrum {
        let mut mz_array = DataArray::from_name_and_type(&ArrayType::MZArray, BinaryDataArrayType::Float64);
        mz_array.extend(mz).unwrap();

        let mut intensity_array =
            DataArray::from_name_and_type(&ArrayType::IntensityArray, BinaryDataArrayType::Float64);
        intensity_array.extend(intensity).unwrap();

        let mut arrays = BinaryArrayMap::new();
        arrays.add(mz_array);
        arrays.add(intensity_array);

        let mut description = SpectrumDescription {
            signal_continuity: mzdata::spectrum::SignalContinuity::Centroid,
            ..Default::default()
        };
        description.precursor.push(Precursor {
            ions: vec![SelectedIon {
                mz: precursor_mz,
                charge: Some(precursor_charge),
                ..Default::default()
            }],
            ..Default::default()
        });

        RawSpectrum { description, arrays }
    }

    /// A profile (unpicked) spectrum must be rejected by [`extract_peaks`], not silently
    /// treated as if every raw data point were its own peak.
    #[test]
    fn test_extract_peaks_rejects_profile_spectrum() {
        let mut spectrum = build_raw_spectrum(&[100.0, 100.1, 100.2], &[1.0, 2.0, 1.0], 2, 300.0);
        spectrum.description.signal_continuity = mzdata::spectrum::SignalContinuity::Profile;

        assert!(super::extract_peaks::<
            mzpeaks::CentroidPeak,
            mzpeaks::DeconvolutedPeak,
            _,
        >(&spectrum)
        .is_err());
    }

    /// Mean squared error between two equal-length score vectors, matching
    /// `ndarray_stats::DeviationExt::mean_sq_err`'s semantics (no square root - despite the
    /// `rmse_*`-style variable names at its call sites, that's what the original 0.0002
    /// thresholds were calibrated against).
    pub fn mean_sq_err(a: &[f64], b: &[f64]) -> f64 {
        assert_eq!(a.len(), b.len());
        let sum_sq: f64 = a.iter().zip(b.iter()).map(|(x, y)| (x - y).powi(2)).sum();
        sum_sq / a.len() as f64
    }

    /// One-time fixture migration: extracts every MS2 spectrum from the source mzML (via
    /// `mzdata`) and writes each one out as `test_files/spectra/scan_<N>.postcard.zst`
    /// (postcard-serialized, then zstd-compressed), independent of whatever fixtures already
    /// exist in that directory. Requires
    /// `LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.mzML` (gitignored, ~1.5GB) at the repo root.
    #[test]
    #[ignore = "Fixture migration, requires the source mzML file locally."]
    fn spectrum_extraction() {
        use mzdata::prelude::*;
        use mzdata::MZReader;

        let reader: MZReader<std::fs::File> =
            MZReader::open_path("LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.mzML").unwrap();

        for spectrum in reader {
            if spectrum.ms_level() != 2 {
                continue;
            }

            // Native IDs look like "controllerType=0 controllerNumber=1 scan=<N>"; `get_spectrum`
            // and the existing fixture naming convention key off the trailing scan number alone.
            let scan = spectrum.id().rsplit("scan=").next().unwrap();

            let bytes = postcard::to_stdvec(&spectrum).unwrap();
            let compressed = zstd::encode_all(bytes.as_slice(), 19).unwrap();
            std::fs::write(
                format!("test_files/spectra/scan_{scan}.postcard.zst"),
                compressed,
            )
            .unwrap();
        }
    }
}
