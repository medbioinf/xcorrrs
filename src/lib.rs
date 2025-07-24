pub mod binning;
pub mod configuration;
pub mod error;
/// Fast xcorr implementation. Less accurate than the reported xcorr in the test data (RSME 0.00018), but faster.
pub mod fast_xcorr;
pub mod scoring_result;
/// Correlation based xcorr. This is closer to the xcorr reported to in the test data (RSME 0.00003).
pub mod xcorr;
// Various utilities
pub mod utils;

/// +/- m/z shift for the xcorr calculation.
pub const BIN_SHIFT: usize = 75;
