pub mod binning;
pub mod configuration;
pub mod error;
#[cfg(feature = "do-not-use-fast-xcorr")]
/// Fast xcorr implementation
pub mod fast_xcorr;
pub mod scoring_result;
/// Correlation based xcorr
pub mod xcorr;
// Various utilities
pub mod utils;
