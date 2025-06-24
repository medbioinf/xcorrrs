use ndarray::ShapeError;
use rustyms::error::CustomError as RustyMsError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum Error {
    #[error("Invalid peptide sequence: {0}")]
    InvalidPeptideSequence(RustyMsError),
    #[error("Invalid charge state: {0}")]
    Padding(ShapeError),
    #[error("Cannot convert {0} to Float when {1}")]
    F64ToFloat(f64, String),
    #[error("Cannot convert {0} to f64 when {1}")]
    FloatToF64(String, String),
    #[error("Cannot convert {0} to usize when {1}")]
    FloatToUsize(String, String),
    #[error("Empty theoretical spectrum")]
    EmptyTheoreticalSpectrum,
    #[error("Empty experimental m/z")]
    EmptyExperimentalSpectrum,
    #[error("m/z ({0}) and intensities ({1}) arrays must have the same length")]
    ExperimentalSpectrumShape(usize, usize),
    #[error("No charge provided via configuration or charge states")]
    NoChargeProvided,
    #[error("Something went wrong: {0}")]
    Correleation(Box<dyn std::error::Error>),
}
