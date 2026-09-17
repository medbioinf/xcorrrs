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
use serde::{Deserialize, Serialize};
use thiserror::Error;

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

    use ndarray::Array1;
    use polars::prelude::*;

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

    /// Reads a spectrum from a Parquet file.
    ///
    /// # Arguments
    /// * `scan` - The scan number to read the spectrum for.
    ///
    pub fn get_spectrum(scan: &str) -> (Array1<f64>, Array1<f64>) {
        let spec_df = ParquetReader::new(
            std::fs::File::open(format!("test_files/spectra/scan_{scan}.parquet")).unwrap(),
        )
        .read_parallel(ParallelStrategy::None)
        .finish()
        .unwrap();

        let mz_array = spec_df["mz"]
            .f64()
            .unwrap()
            .to_ndarray()
            .unwrap()
            .to_owned();

        let intensity_array = spec_df["intensity"]
            .f64()
            .unwrap()
            .to_ndarray()
            .unwrap()
            .to_owned();

        (mz_array, intensity_array)
    }

    // TODO: Replace with omicstools with mzdata and make it independent of comet's test data
    // /// Extracts spectru, data arrays from the mzML file and saves them as Parquet files without any further metadata
    // /// Important do not set TEST_NUMBER_OF_PSMS environment variable, as this will limit the spectrum extraction to a certain number of scans.
    // ///
    // #[test]
    // #[ignore = "Spectrum extration."]
    // fn spectrum_extraction() {
    //     let comet_df = read_test_data();

    //     let mut mzml_byte_reader = std::io::BufReader::new(
    //         std::fs::File::open("test_files/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.mzML")
    //             .unwrap(),
    //     );
    //     let mut mzml = dihardts_omicstools::proteomics::io::mzml::reader::Reader::read_indexed(
    //         &mut mzml_byte_reader,
    //         None,
    //         true,
    //         false,
    //     )
    //     .unwrap();

    //     for i in 0..comet_df.height() {
    //         let scan = comet_df["scan"].i64().unwrap().get(i).unwrap();

    //         let binary_data_array_list = mzml
    //             .get_spectrum(&format!("controllerType=0 controllerNumber=1 scan={scan}"))
    //             .unwrap()
    //             .binary_data_array_list;

    //         let mz_array = binary_data_array_list
    //             .get_mz_array()
    //             .unwrap()
    //             .deflate_data()
    //             .unwrap();
    //         let intensity_array = binary_data_array_list
    //             .get_intensity_array()
    //             .unwrap()
    //             .deflate_data()
    //             .unwrap();

    //         let mut spec_frame = DataFrame::new(vec![
    //             Column::new("mz".into(), mz_array),
    //             Column::new("intensity".into(), intensity_array),
    //         ])
    //         .unwrap();

    //         let writer = ParquetWriter::new(
    //             std::fs::File::create(format!("test_files/spectra/scan_{scan}.parquet")).unwrap(),
    //         )
    //         .with_compression(ParquetCompression::Zstd(Some(
    //             ZstdLevel::try_new(22).unwrap(),
    //         )));
    //         writer.finish(&mut spec_frame).unwrap();
    //     }
    // }
}
