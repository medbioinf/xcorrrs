use context_error::{BasicKind, BoxedError};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum Error {
    #[error("Empty experimental m/z")]
    EmptyExperimentalSpectrum,
    #[error("m/z ({0}) and intensities ({1}) arrays must have the same length")]
    ExperimentalSpectrumShape(usize, usize),
    #[error("Invalid configuration: {0}")]
    InvalidConfiguration(String),
    #[error("Invalid ProForma string: {0:?}")]
    InvalidProForma(Vec<BoxedError<'static, BasicKind>>),
    #[error("Unable to create theoretical spectra: {0:?}")]
    TheoreticalSpectrum(Box<crate::utils::Error>),
}
