use context_error::{BasicKind, BoxedError};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum Error {
    #[error("Empty experimental m/z")]
    EmptyExperimentalSpectrum,
    #[error("Spectrum has no precursor charge")]
    MissingCharge,
    #[error(transparent)]
    InvalidPrecursor(#[from] crate::utils::PrecursorMultiplicityError),
    #[error("Spectrum is not centroided (profile signal continuity) - pick peaks before scoring")]
    UncentroidedSpectrum,
    #[error("Invalid configuration: {0}")]
    InvalidConfiguration(String),
    #[error("Invalid ProForma string: {0:?}")]
    InvalidProForma(Vec<BoxedError<'static, BasicKind>>),
    #[error("Unable to create theoretical spectra: {0:?}")]
    TheoreticalSpectrum(Box<crate::utils::Error>),
}
