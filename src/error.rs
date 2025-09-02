use rustyms::error::CustomError as RustyMsError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum Error {
    #[error("Invalid peptide sequence: {0}")]
    InvalidPeptideSequence(RustyMsError),
    #[error("Empty experimental m/z")]
    EmptyExperimentalSpectrum,
    #[error("m/z ({0}) and intensities ({1}) arrays must have the same length")]
    ExperimentalSpectrumShape(usize, usize),
    #[error("Invalid configuration: {0}")]
    InvalidConfiguration(String),
}
