use std::collections::{HashMap, HashSet};

use context_error::{BasicKind, BoxedError};
use mzannotate::{annotation::model::FragmentationModel, fragment::FragmentType};
use mzcore::{
    chemistry::{AmbiguousMolecule, MassMode, Molecule},
    molecular_formula,
    ontology::Ontologies,
    quantities::Tolerance,
    sequence::{IsAminoAcid, Modification, PeptidoformIon},
    system::{mass_over_charge::thomson, MassOverCharge},
};
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::{
    msamanda::deconvolution::{get_mass_index, Spectrum},
    utils::{build_fragmentation_model, create_theoretical_fragments, Ion},
};

/// Configuration for the MS Amanda scoring algorithm.
///
/// See: <https://pubs.acs.org/doi/10.1021/pr500202e>
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Configuration {
    /// Fragment ion mass tolerance, used both to match theoretical ions against picked peaks
    /// and to strip the precursor peak from the spectrum.
    pub fragment_tolerance: Tolerance<MassOverCharge>,
    /// Maximum charge state to consider for theoretical fragment ions.
    pub max_fragment_charge: usize,
    /// Ions to use in the fragmentation model.
    pub ions: Vec<Ion>,
    /// Whether to include water/ammonia neutral-loss variants of b/y ions.
    pub use_neutral_loss_ions: bool,
    /// Whether to consider the methionine-sulfoxide neutral loss (-CH4OS, ~64 Da) as an extra
    /// candidate ion for any b/y fragment that spans an oxidized methionine residue. This is a
    /// collision-induced side-chain loss (vibrational excitation), so it is only chemically
    /// meaningful for CID/HCD backbone fragmentation (b/y ions) and is never applied to c/z
    /// (ETD/ECD) ions regardless of this setting. Off by default; enable explicitly for CID/HCD
    /// searches where oxidized-Met PSMs matter.
    pub use_methionine_sulfoxide_neutral_loss: bool,
    /// Smallest peak-picking depth (peaks kept per window) to try.
    pub min_m: usize,
    /// Largest peak-picking depth to try. MS Amanda sweeps `min_m..=max_m` and keeps whichever
    /// depth yields the highest score.
    pub max_m: usize,
    /// Width, in Da, of each peak-picking window.
    pub peak_picking_window: f64,
    /// Whether to strip peaks within `fragment_tolerance` of the precursor m/z before scoring.
    pub remove_precursor: bool,
}

impl Default for Configuration {
    /// Default is for high mass accuracy MS2 data, matching the paper's recommended setup.
    fn default() -> Self {
        Self {
            fragment_tolerance: Tolerance::new_ppm(10.0),
            max_fragment_charge: 3,
            ions: vec![Ion::B, Ion::Y],
            use_neutral_loss_ions: false,
            use_methionine_sulfoxide_neutral_loss: false,
            min_m: 1,
            max_m: 10,
            peak_picking_window: 100.0,
            remove_precursor: true,
        }
    }
}

/// Finalized configuration for the MS Amanda algorithm, carrying the built `FragmentationModel`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FinalizedConfiguration {
    /// Fragmentation model built from the selected ions.
    pub fragmentation_model: FragmentationModel,
    /// Original configuration that was used to create this finalized configuration.
    configuration: Configuration,
}

impl std::ops::Deref for FinalizedConfiguration {
    type Target = Configuration;

    fn deref(&self) -> &Self::Target {
        &self.configuration
    }
}

impl From<Configuration> for FinalizedConfiguration {
    fn from(configuration: Configuration) -> Self {
        let fragmentation_model =
            build_fragmentation_model(&configuration.ions, configuration.use_neutral_loss_ions);

        Self {
            fragmentation_model,
            configuration,
        }
    }
}

impl Default for FinalizedConfiguration {
    fn default() -> Self {
        Configuration::default().into()
    }
}

#[derive(Error, Debug)]
pub enum Error {
    #[error("Empty spectrum")]
    EmptySpectrum,
    #[error("Spectrum has no precursor charge and no charge fallback was given")]
    MissingCharge,
    #[error("Spectrum has no precursor m/z")]
    MissingPrecursorMz,
    #[error("Spectrum is not centroided (profile signal continuity) - pick peaks before scoring")]
    UncentroidedSpectrum,
    #[error("Invalid ProForma string: {0:?}")]
    InvalidProForma(Vec<BoxedError<'static, BasicKind>>),
    #[error("Unable to create theoretical spectra: {0:?}")]
    TheoreticalSpectrum(Box<crate::utils::Error>),
}

/// Result of scoring a single peptide against a spectrum with MS Amanda.
#[derive(Debug, Clone, PartialEq)]
pub struct MSAmandaScoringResult {
    /// Eq (7): `-10 * log10(p_value / explained_ion_current)`.
    pub score: f64,
    /// Eq (1): cumulative binomial probability of matching `ions_matched` picked peaks by chance.
    pub p_value: f64,
    /// Eq (6): fraction of picked-peak intensity explained by matched peaks (`eif`).
    pub explained_ion_current: f64,
    /// Number of picked peaks (`N`) at the winning peak-picking depth.
    pub picked_peaks: usize,
    /// Winning peak-picking depth (`m`), out of the configured sweep range.
    pub peak_picking_depth: usize,
    /// Number of distinct picked peaks matched by at least one theoretical ion (`n`), the
    /// quantity that feeds into the binomial probability and `eif`.
    pub matched_peaks: usize,
    /// Total number of theoretical fragment ions generated.
    pub ions_total: usize,
    /// Number of theoretical fragment ions that matched at least one picked peak.
    pub ions_matched: usize,
    pub min_theoretical_mass: f64,
    pub max_theoretical_mass: f64,
}

impl MSAmandaScoringResult {
    pub fn round_score(&self, decimals: u16) -> f64 {
        let factor = 10f64.powi(decimals as i32);
        (self.score * factor).round() / factor
    }
}

impl std::fmt::Display for MSAmandaScoringResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Score: {:.2}, p-value: {:.3e}, eif: {:.3}, m: {}, Picked Peaks: {}, Matched Peaks: {}, Ions Total: {}, Ions Matched: {}, Min Mass: {:.2}, Max Mass: {:.2}",
            self.score,
            self.p_value,
            self.explained_ion_current,
            self.peak_picking_depth,
            self.picked_peaks,
            self.matched_peaks,
            self.ions_total,
            self.ions_matched,
            self.min_theoretical_mass,
            self.max_theoretical_mass
        )
    }
}

/// A single picked peak (post peak-picking), decoupled from the spectrum's internal encoding.
#[derive(Debug, Clone, Copy)]
struct PickedPeak {
    mz: f64,
    intensity: f64,
}

/// Spectrum peaks bucketed into peak-picking windows, each sorted by intensity descending, so
/// that peak-picking at any depth `m` is just taking a prefix of each window.
struct WindowedPeaks {
    windows: Vec<Vec<PickedPeak>>,
}

/// Peak picking (paper, preprocessing step before eq 1): buckets every peak of `spectrum` into
/// `window_width`-wide windows (anchored at the spectrum's lowest m/z, 100 Da per paper),
/// sorting each window by intensity descending so [`pick_peaks`] can keep the top `m`.
fn window_peaks(spectrum: &Spectrum, window_width: f64) -> Option<WindowedPeaks> {
    let min_mz = spectrum.masses.iter().next().copied()?.into_inner();

    let mut windows: HashMap<usize, Vec<PickedPeak>> = HashMap::new();

    for &mz in &spectrum.masses {
        let mz_val = mz.into_inner();
        let intensity = spectrum
            .peaks
            .get(&get_mass_index(mz))
            .copied()
            .unwrap_or(OrderedFloat(0.0))
            .into_inner();
        let window_index = ((mz_val - min_mz) / window_width).floor() as usize;
        windows.entry(window_index).or_default().push(PickedPeak {
            mz: mz_val,
            intensity,
        });
    }

    for peaks in windows.values_mut() {
        peaks.sort_by(|a, b| {
            b.intensity
                .partial_cmp(&a.intensity)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    Some(WindowedPeaks {
        windows: windows.into_values().collect(),
    })
}

/// Peak picking (paper, preprocessing step before eq 1): picks the `m` most intense peaks per
/// window, merged and sorted by m/z ascending.
fn pick_peaks(windowed: &WindowedPeaks, m: usize) -> Vec<PickedPeak> {
    let mut picked: Vec<PickedPeak> = windowed
        .windows
        .iter()
        .flat_map(|peaks| peaks.iter().take(m).copied())
        .collect();
    picked.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(std::cmp::Ordering::Equal));
    picked
}

/// Eq (3) numerator (`|F|*2t - O(F)`, using eq (4)-(5)'s `O(F)`): total width covered by the
/// union of (possibly overlapping) `[lower, upper]` theoretical-ion tolerance intervals.
/// Computed directly as an interval-union merge instead of `|F|*2t` minus the paper's sequential
/// pairwise-sum overlap `O(F)`, but gives the same total for sorted, possibly-chained overlaps.
fn union_width(mut intervals: Vec<(f64, f64)>) -> f64 {
    if intervals.is_empty() {
        return 0.0;
    }
    intervals.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut total = 0.0;
    let (mut current_low, mut current_high) = intervals[0];
    for &(low, high) in &intervals[1..] {
        if low > current_high {
            total += current_high - current_low;
            current_low = low;
            current_high = high;
        } else {
            current_high = current_high.max(high);
        }
    }
    total += current_high - current_low;
    total
}

/// Eq (1): `P(n,p,N)`, cumulative binomial probability of matching at least `n` out of
/// `capital_n` trials, each with success probability `p`:
/// `P(X >= n) = sum_{i=n}^{capital_n} C(capital_n,i) p^i (1-p)^(capital_n-i)`.
///
/// Computed via a forward recurrence in log-space (`pmf(i+1) = pmf(i) * (N-i)/(i+1) * p/(1-p)`)
/// instead of evaluating binomial coefficients directly, to stay numerically stable for large
/// `capital_n`. Crucially, the upper tail is summed *directly* (not as `1 - P(X < n)`) since for
/// strong matches `P(X < n)` rounds to `1.0` in `f64`, which would make the complement collapse
/// to exactly `0.0` and destroy the (very informative) tiny p-value.
fn cumulative_binomial_at_least(n: usize, capital_n: usize, p: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    if n > capital_n {
        return 0.0;
    }
    if p <= 0.0 {
        return 0.0;
    }
    if p >= 1.0 {
        return 1.0;
    }

    // log(pmf(0)) = N * ln(1-p); advance the recurrence up to i = n without accumulating,
    // since only the i >= n terms are needed.
    let mut log_pmf = (capital_n as f64) * (1.0 - p).ln();
    for i in 0..n {
        let remaining = (capital_n - i) as f64;
        let taken = (i + 1) as f64;
        log_pmf += remaining.ln() - taken.ln() + p.ln() - (1.0 - p).ln();
    }

    let mut cdf_at_least_n = 0.0;
    for i in n..=capital_n {
        cdf_at_least_n += log_pmf.exp();
        if i < capital_n {
            let remaining = (capital_n - i) as f64;
            let taken = (i + 1) as f64;
            log_pmf += remaining.ln() - taken.ln() + p.ln() - (1.0 - p).ln();
        }
    }

    cdf_at_least_n.clamp(0.0, 1.0)
}

/// MS Amanda scorer for a single (deisotoped) experimental spectrum.
pub struct MSAmanda<'a> {
    config: &'a FinalizedConfiguration,
    spectrum: Spectrum,
    fragment_charge: usize,
}

impl MSAmanda<'_> {
    /// Creates a new MS Amanda scorer for `spectrum`, which is expected to already be
    /// deisotoped (see [`crate::msamanda::deconvolution::deisotope_spectrum`]). Precursor peaks
    /// are stripped here if `config.remove_precursor` is set.
    pub fn new<'a>(
        config: &'a FinalizedConfiguration,
        mut spectrum: Spectrum,
    ) -> Result<MSAmanda<'a>, Error> {
        if spectrum.masses.is_empty() {
            return Err(Error::EmptySpectrum);
        }

        if config.remove_precursor {
            Self::remove_precursor_peaks(&mut spectrum, config.fragment_tolerance);
        }

        if spectrum.masses.is_empty() {
            return Err(Error::EmptySpectrum);
        }

        let mut fragment_charge = (spectrum.precursor.charge - 1).max(1) as usize;
        if fragment_charge > config.max_fragment_charge {
            fragment_charge = config.max_fragment_charge;
        }

        Ok(MSAmanda {
            config,
            spectrum,
            fragment_charge,
        })
    }

    /// Removes every peak within `tolerance` of the spectrum's precursor m/z.
    fn remove_precursor_peaks(spectrum: &mut Spectrum, tolerance: Tolerance<MassOverCharge>) {
        let precursor_mz = MassOverCharge::new::<thomson>(spectrum.precursor.mz.into_inner());
        let (lower, upper) = tolerance.bounds(precursor_mz);
        let (lower, upper) = (lower.value, upper.value);

        let to_remove: Vec<OrderedFloat<f64>> = spectrum
            .masses
            .iter()
            .copied()
            .filter(|mz| {
                let value = mz.into_inner();
                value >= lower && value <= upper
            })
            .collect();

        for mz in to_remove {
            spectrum.masses.remove(&mz);
            spectrum.peaks.remove(&get_mass_index(mz));
        }
    }

    /// Scores `peptide` against the spectrum this scorer was built for.
    ///
    /// Sweeps peak-picking depth `m` over `config.min_m..=config.max_m` and returns the
    /// highest-scoring result, per the paper's algorithm.
    pub fn score_peptide(
        &self,
        peptide: &str,
        ontologies: Option<&Ontologies>,
    ) -> Result<MSAmandaScoringResult, Error> {
        let empty_ontology = Ontologies::empty();
        let ontologies = ontologies.unwrap_or(&empty_ontology);
        let (peptidoform_ion, parsing_errors) = PeptidoformIon::pro_forma(peptide, ontologies)
            .map_err(|err| {
                Error::InvalidProForma(err.into_iter().map(BoxedError::to_owned).collect())
            })?;

        if !parsing_errors.is_empty() {
            return Err(Error::InvalidProForma(
                parsing_errors
                    .into_iter()
                    .map(BoxedError::to_owned)
                    .collect(),
            ));
        }

        let (min_theoretical_mass, max_theoretical_mass) =
            match peptidoform_ion.formulas().mass_bounds().into_option() {
                Some((min, max)) => (min.monoisotopic_mass().value, max.monoisotopic_mass().value),
                None => (-1.0, -1.0),
            };

        let theoretical_fragments = create_theoretical_fragments(
            &peptidoform_ion,
            &self.config.fragmentation_model,
            self.fragment_charge,
        )
        .map_err(|err| Error::TheoreticalSpectrum(Box::new(err)))?;

        // 0-based residue indices carrying an oxidation (+15.9949 Da) on methionine, however that
        // modification happens to be encoded (a bare mass delta, a `Database`/UNIMOD:35 lookup,
        // a molecular formula, ...) - `Molecule::monoisotopic_mass` resolves all of them to the
        // same underlying mass, so detection isn't tied to one particular ProForma spelling.
        // Real MS Amanda runs (per `settings.xml`'s `MaxNumberNeutralLossModifications`) consider
        // the characteristic methionine-sulfoxide neutral loss (-CH4OS, ~64 Da) as an extra
        // candidate ion for any b/y fragment spanning such a residue, on top of the unmodified ion.
        let oxidized_met_positions: Vec<usize> = if self.config.use_methionine_sulfoxide_neutral_loss
        {
            peptidoform_ion
                .singular_ref()
                .map(|peptidoform| {
                    peptidoform
                        .sequence()
                        .iter()
                        .enumerate()
                        .filter(|(_, element)| {
                            element.aminoacid.one_letter_code() == Some('M')
                                && element.modifications.iter().any(|modification| {
                                    let Modification::Simple(simple) = modification else {
                                        return false;
                                    };
                                    (simple.monoisotopic_mass().mass().value - 15.994_915).abs()
                                        < 0.01
                                })
                        })
                        .map(|(index, _)| index)
                        .collect()
                })
                .unwrap_or_default()
        } else {
            Vec::new()
        };
        let met_sulfoxide_loss_mass = molecular_formula!(C 1 H 4 O 1 S 1).monoisotopic_mass().value;

        // Each theoretical ion's raw m/z plus its match-tolerance window; both are constant
        // across the `m` sweep below, so are computed once up front. The raw m/z is kept
        // alongside the window so that, per depth `m`, the ions can be restricted to `F(s,pep,m)`
        // (see `Configuration` docs / MS Amanda paper eq. 2): only ions within the tolerance-
        // expanded range of the picked peaks contribute to the chance-match probability `p`.
        let mut theoretical_ions: Vec<(f64, f64, f64)> = Vec::with_capacity(
            theoretical_fragments.len() + oxidized_met_positions.len(),
        );
        for fragment in &theoretical_fragments {
            let Some(mz) = fragment.mz(MassMode::Monoisotopic) else {
                continue;
            };
            let (lower, upper) = self.config.fragment_tolerance.bounds(mz);
            theoretical_ions.push((mz.value, lower.value, upper.value));

            let covers_oxidized_met = match &fragment.ion {
                FragmentType::b(position, _) => oxidized_met_positions
                    .iter()
                    .any(|&index| index < position.series_number),
                FragmentType::y(position, _) => oxidized_met_positions.iter().any(|&index| {
                    index >= position.sequence_length.saturating_sub(position.series_number)
                }),
                _ => false,
            };

            if covers_oxidized_met {
                let loss_mz = met_sulfoxide_loss_mass / fragment.charge.value as f64;
                let loss_mz_value = mz.value - loss_mz;
                if loss_mz_value > 0.0 {
                    let loss_mz = MassOverCharge::new::<thomson>(loss_mz_value);
                    let (lower, upper) = self.config.fragment_tolerance.bounds(loss_mz);
                    theoretical_ions.push((loss_mz_value, lower.value, upper.value));
                }
            }
        }

        let ions_total = theoretical_ions.len();

        let windowed = window_peaks(&self.spectrum, self.config.peak_picking_window)
            .ok_or(Error::EmptySpectrum)?;

        let mut best: Option<MSAmandaScoringResult> = None;

        // Eq (7): sweep peak-picking depth m over [1..10] (config min_m..=max_m), scoring each
        // depth independently and keeping the max-scoring one below.
        for m in self.config.min_m..=self.config.max_m {
            let picked = pick_peaks(&windowed, m);
            let capital_n = picked.len();
            if capital_n == 0 {
                continue;
            }

            // e1(s,m)-t and eN(s,m)+t (eq 2/3 denominator): tolerance-expanded m/z range of the
            // first and last picked peak.
            let (e1_lower, _) = self
                .config
                .fragment_tolerance
                .bounds(MassOverCharge::new::<thomson>(picked.first().unwrap().mz));
            let (_, en_upper) = self
                .config
                .fragment_tolerance
                .bounds(MassOverCharge::new::<thomson>(picked.last().unwrap().mz));
            let e1_lower = e1_lower.value;
            let en_upper = en_upper.value;
            let range = en_upper - e1_lower;

            // Eq (2): F(s,pep,m), theoretical ions whose raw m/z falls within [e1_lower, en_upper];
            // union_width of those ions' tolerance windows gives the eq (3) numerator directly.
            let overlap_width = union_width(
                theoretical_ions
                    .iter()
                    .filter(|&&(raw_mz, _, _)| raw_mz >= e1_lower && raw_mz <= en_upper)
                    .map(|&(_, lower, upper)| (lower, upper))
                    .collect(),
            );

            let total_intensity: f64 = picked.iter().map(|peak| peak.intensity).sum();

            let mut matched_indices: HashSet<usize> = HashSet::new();
            let mut ions_matched = 0usize;

            for &(_, lower, upper) in &theoretical_ions {
                let mut matched_this_ion = false;
                for (idx, peak) in picked.iter().enumerate() {
                    if peak.mz >= lower && peak.mz <= upper {
                        matched_indices.insert(idx);
                        matched_this_ion = true;
                    }
                }
                if matched_this_ion {
                    ions_matched += 1;
                }
            }

            // n, matched peak count, feeds eq (1) and eq (6).
            let n = matched_indices.len();
            let matched_intensity: f64 = matched_indices
                .iter()
                .map(|&idx| picked[idx].intensity)
                .sum();
            // Eq (6): eif(s,pep,m) = sum I(matched) / sum I(picked).
            let eif = if total_intensity > 0.0 {
                matched_intensity / total_intensity
            } else {
                0.0
            };

            // Eq (3): p(s,pep,m) = (|F|*2t - O(F)) / ((eN+t)-(e1-t)) = overlap_width / range.
            let p = if n == 0 || range <= 0.0 {
                0.0
            } else {
                (overlap_width / range).clamp(0.0, 1.0)
            };

            // Eq (1): P(n,p,N) with N = capital_n (picked peaks), n = matched peaks, p from eq (3).
            let p_value = cumulative_binomial_at_least(n, capital_n, p);

            // Eq (7): S(s,pep) = -10 * log10(P(n,p,N) / eif) for this depth m.
            let score = if n == 0 || eif <= 0.0 {
                0.0
            } else {
                -10.0 * (p_value.max(1e-300) / eif).log10()
            };

            let candidate = MSAmandaScoringResult {
                score,
                p_value,
                explained_ion_current: eif,
                picked_peaks: capital_n,
                peak_picking_depth: m,
                matched_peaks: n,
                ions_total,
                ions_matched,
                min_theoretical_mass,
                max_theoretical_mass,
            };

            if best
                .as_ref()
                .is_none_or(|current_best| candidate.score > current_best.score)
            {
                best = Some(candidate);
            }
        }

        best.ok_or(Error::EmptySpectrum)
    }
}

#[cfg(test)]
mod tests {
    use ordered_float::OrderedFloat;
    use std::{
        collections::{BTreeSet, HashMap},
        path::PathBuf,
    };

    use super::*;
    use crate::msamanda::deconvolution::{self, Precursor};

    /// Adapts `Configuration::fragment_tolerance` (an `mzcore` ppm/absolute tolerance) to the
    /// local [`deconvolution::Tolerance`] trait expected by [`deconvolution::deisotope_spectrum`].
    struct FragmentToleranceAdapter<'a>(&'a Tolerance<MassOverCharge>);

    impl deconvolution::Tolerance for FragmentToleranceAdapter<'_> {
        fn calculate_lower_bound(&self, mass: OrderedFloat<f64>) -> OrderedFloat<f64> {
            let (lower, _) = self.0.bounds(MassOverCharge::new::<thomson>(mass.into_inner()));
            OrderedFloat(lower.value)
        }

        fn calculate_upper_bound(&self, mass: OrderedFloat<f64>) -> OrderedFloat<f64> {
            let (_, upper) = self.0.bounds(MassOverCharge::new::<thomson>(mass.into_inner()));
            OrderedFloat(upper.value)
        }
    }

    fn build_spectrum(peaks: &[(f64, f64)], precursor_charge: i32, precursor_mz: f64) -> Spectrum {
        let mut peak_map = HashMap::new();
        let mut masses = BTreeSet::new();

        for &(mz, intensity) in peaks {
            let mz = OrderedFloat(mz);
            masses.insert(mz);
            peak_map.insert(get_mass_index(mz), OrderedFloat(intensity));
        }

        Spectrum {
            peaks: peak_map,
            masses,
            precursor: Precursor {
                charge: precursor_charge,
                mz: OrderedFloat(precursor_mz),
            },
        }
    }

    #[test]
    fn test_cumulative_binomial_edge_cases() {
        assert_eq!(cumulative_binomial_at_least(0, 10, 0.1), 1.0);
        assert_eq!(cumulative_binomial_at_least(11, 10, 0.1), 0.0);
    }

    #[test]
    fn test_cumulative_binomial_known_values() {
        // N=2, p=0.5: P(X>=1) = 0.75, P(X>=2) = 0.25 (hand-computed from the binomial pmf).
        let p_ge_1 = cumulative_binomial_at_least(1, 2, 0.5);
        let p_ge_2 = cumulative_binomial_at_least(2, 2, 0.5);
        assert!((p_ge_1 - 0.75).abs() < 1e-9);
        assert!((p_ge_2 - 0.25).abs() < 1e-9);
    }

    #[test]
    fn test_peak_picking_depth() {
        // Two peaks in the same 100 Da window; depth 1 keeps only the most intense.
        let spectrum = build_spectrum(&[(100.0, 5.0), (110.0, 50.0)], 2, 300.0);
        let windowed = window_peaks(&spectrum, 100.0).unwrap();

        let picked_1 = pick_peaks(&windowed, 1);
        assert_eq!(picked_1.len(), 1);
        assert_eq!(picked_1[0].mz, 110.0);

        let picked_2 = pick_peaks(&windowed, 2);
        assert_eq!(picked_2.len(), 2);
    }

    #[test]
    fn test_precursor_peak_removed() {
        let mut spectrum = build_spectrum(&[(100.0, 5.0), (300.0, 999.0)], 2, 300.0);
        let config: FinalizedConfiguration = Configuration {
            fragment_tolerance: Tolerance::new_absolute(MassOverCharge::new::<thomson>(0.5)),
            ..Configuration::default()
        }
        .into();

        MSAmanda::remove_precursor_peaks(&mut spectrum, config.fragment_tolerance);

        assert_eq!(spectrum.masses.len(), 1);
        assert!(spectrum.masses.contains(&OrderedFloat(100.0)));
    }

    /// Builds ProForma string. for the test data we use oxidation of M and carbamidomethylation of C.
    /// Need replacement of m by M[+15.9949] ad C[+57.02146]
    fn build_proforma(sequence: &str) -> String {
        let mut proforma_string = sequence.replace("m", "M[+15.994915]").replace("c", "C");
        if proforma_string.contains("C") {
            proforma_string = format!("<[+57.021464]@C>{proforma_string}")
        }
        proforma_string
    }

    /// Reads MS Amanda's own search results for the same spectra used by `xcorr`'s
    /// `test_xcorr` (see `crate::utils::tests::read_test_data`), for comparison against this
    /// implementation. `TEST_NUMBER_OF_PSMS` limits to the X best scored PSMs, same convention
    /// as the Comet reader.
    fn read_test_data() -> polars::prelude::DataFrame {
        use polars::prelude::*;

        let mut msamanda_df = CsvReadOptions::default()
            .with_has_header(true)
            .with_parse_options(
                CsvParseOptions::default()
                    .with_separator(b'\t')
                    .with_comment_prefix(Some("#")),
            )
            .try_into_reader_with_file_path(Some(PathBuf::from(
                "test_files/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.msamanda.tsv",
            )))
            .unwrap()
            .finish()
            .unwrap();

        msamanda_df
            .sort_in_place(
                ["Amanda Score"],
                SortMultipleOptions::default().with_order_descending(true),
            )
            .unwrap();

        match std::env::var("TEST_NUMBER_OF_PSMS") {
            Ok(number_of_psms) => msamanda_df.slice(0, number_of_psms.parse::<usize>().unwrap()),
            Err(_) => msamanda_df,
        }
    }

    // Test the MS Amanda scoring implementation against real MS Amanda 3.0 output, using the
    // same spectra as `xcorr::fast_xcorr::tests::test_xcorr`.
    #[test]
    fn test_msamanda_scoring() {
        use polars::prelude::*;
        use rayon::prelude::*;

        use crate::utils::tests::{get_spectrum, mean_sq_err};

        let msamanda_df = read_test_data();

        let config: FinalizedConfiguration = Configuration {
            // Matches `test_files/settings.xml`'s `<MS2Tol Unit="Da">0.02</MS2Tol>`, the real
            // MS Amanda run's fragment tolerance used to generate the ground-truth TSV.
            fragment_tolerance: Tolerance::new_absolute(MassOverCharge::new::<thomson>(0.02)),
            // Matches `test_files/settings.xml`'s `<MaxAllowedChargeState>4+</MaxAllowedChargeState>`.
            max_fragment_charge: 4,
            // This dataset is a CID/HCD (b/y ion) search per `settings.xml`'s
            // `<Instrument>b, y</Instrument>`, and its `MaxNumberNeutralLossModifications`
            // setting means modification-specific losses like methionine-sulfoxide were
            // considered by the real run.
            use_methionine_sulfoxide_neutral_loss: true,
            ..Configuration::default()
        }
        .into();

        #[allow(clippy::type_complexity)]
        let results: Vec<(i64, String, f64, f64, f64, f64, f64, f64, u64, u64)> = (0..msamanda_df
            .height())
            .into_par_iter()
            .filter_map(|idx| {
                let scan = msamanda_df["Scan Number"].i64().unwrap().get(idx).unwrap();

                // Not every scan in the MS Amanda results has a pre-extracted spectrum fixture
                // file (see `crate::utils::tests::spectrum_extraction`); skip those.
                if !PathBuf::from(format!("test_files/spectra/scan_{scan}.postcard.zst")).exists()
                {
                    return None;
                }

                let amanda_score = msamanda_df["Amanda Score"].f64().unwrap().get(idx).unwrap();
                let binom_score = msamanda_df["binom score"].f64().unwrap().get(idx).unwrap();
                let assigned_intensity_fraction = msamanda_df["assigned intensity fraction"]
                    .f64()
                    .unwrap()
                    .get(idx)
                    .unwrap();
                let comet_matched_peaks = msamanda_df["Nr of matched peaks"]
                    .i64()
                    .unwrap()
                    .get(idx)
                    .unwrap() as u64;
                let sequence = msamanda_df["Sequence"].str().unwrap().get(idx).unwrap();
                let charge = msamanda_df["Charge"].i64().unwrap().get(idx).unwrap() as i32;

                let proforma_peptide = build_proforma(sequence);
                let spectrum = get_spectrum(scan.to_string().as_str());
                let mut spectrum =
                    deconvolution::Spectrum::from_mzdata_spectrum(&spectrum, Some(charge as usize))
                        .ok()?;
                deconvolution::deisotope_spectrum(
                    &mut spectrum,
                    &FragmentToleranceAdapter(&config.fragment_tolerance),
                );

                let scoring = MSAmanda::new(&config, spectrum).ok()?;
                let result = scoring.score_peptide(&proforma_peptide, None).ok()?;

                Some((
                    scan,
                    proforma_peptide,
                    amanda_score,
                    result.round_score(3),
                    binom_score,
                    result.p_value,
                    assigned_intensity_fraction,
                    result.explained_ion_current,
                    comet_matched_peaks,
                    result.matched_peaks as u64,
                ))
            })
            .collect();

        assert!(
            !results.is_empty(),
            "no MS Amanda PSMs could be scored - check test_files/spectra availability"
        );

        let scan_col: Vec<i64> = results.iter().map(|r| r.0).collect();
        let peptide_col: Vec<String> = results.iter().map(|r| r.1.clone()).collect();
        let amanda_score_col: Vec<f64> = results.iter().map(|r| r.2).collect();
        let rs_score_col: Vec<f64> = results.iter().map(|r| r.3).collect();
        let amanda_binom_col: Vec<f64> = results.iter().map(|r| r.4).collect();
        let rs_p_value_col: Vec<f64> = results.iter().map(|r| r.5).collect();
        let amanda_eif_col: Vec<f64> = results.iter().map(|r| r.6).collect();
        let rs_eif_col: Vec<f64> = results.iter().map(|r| r.7).collect();
        let amanda_matched_col: Vec<u64> = results.iter().map(|r| r.8).collect();
        let rs_matched_col: Vec<u64> = results.iter().map(|r| r.9).collect();

        let mut comparison_df = DataFrame::new(vec![
            Column::new("scan".into(), scan_col.clone()),
            Column::new("proforma_peptide".into(), peptide_col),
            Column::new("amanda_score".into(), amanda_score_col.clone()),
            Column::new("rs_score".into(), rs_score_col.clone()),
            Column::new("amanda_binom_score".into(), amanda_binom_col),
            Column::new("rs_p_value".into(), rs_p_value_col),
            Column::new("amanda_eif".into(), amanda_eif_col),
            Column::new("rs_eif".into(), rs_eif_col),
            Column::new("amanda_matched_peaks".into(), amanda_matched_col),
            Column::new("rs_matched_peaks".into(), rs_matched_col),
        ])
        .unwrap();

        CsvWriter::new(std::fs::File::create("msamanda_comparison.tsv").unwrap())
            .with_separator(b'\t')
            .finish(&mut comparison_df)
            .unwrap();

        // Normalize both score columns by their max, like `xcorr::fast_xcorr::tests::test_xcorr`
        // does, since the absolute scale of a from-the-paper reimplementation need not match the
        // original tool exactly.
        let max_amanda_score = amanda_score_col
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let max_rs_score = rs_score_col
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let max_score = max_amanda_score.max(max_rs_score);

        if std::env::var("VERBOSE").is_ok() {
            let mut plot = plotly::Plot::new();
            let diagonal_trace =
                plotly::Scatter::new(vec![0.0, max_amanda_score], vec![0.0, max_rs_score])
                    .mode(plotly::common::Mode::Lines)
                    .marker(plotly::common::Marker::default().color("red"))
                    .hover_info(plotly::common::HoverInfo::None)
                    .show_legend(false);

            let correlation_trace =
                plotly::Scatter::new(amanda_score_col.clone(), rs_score_col.clone())
                    .hover_text_array(
                        scan_col
                            .iter()
                            .map(|scan| format!("scan = {}", scan))
                            .collect::<Vec<String>>(),
                    )
                    .mode(plotly::common::Mode::Markers)
                    .marker(plotly::common::Marker::default().color("blue"))
                    .show_legend(false);

            plot.add_trace(diagonal_trace);
            plot.add_trace(correlation_trace);

            plot.set_layout(
                plotly::Layout::new()
                    .title("MSAmanda vs mzscoring msamanda score")
                    .x_axis(
                        plotly::layout::Axis::new()
                            .title("MSAmanda score")
                            .constrain(plotly::layout::AxisConstrain::Domain),
                    )
                    .y_axis(
                        plotly::layout::Axis::new()
                            .title("mzscoring msamanda score")
                            .scale_anchor("x"),
                    ),
            );
            plot.write_html("99-msamanda_score_vs_amanda_score.html");
        }

        let scaled_amanda_score = amanda_score_col
            .into_iter()
            .map(|value| value / max_score)
            .collect::<Vec<f64>>();
        let scaled_rs_score = rs_score_col
            .into_iter()
            .map(|value| value / max_score)
            .collect::<Vec<f64>>();

        let rmse = mean_sq_err(&scaled_amanda_score, &scaled_rs_score);

        println!("RMSE amanda_score vs rs_score (normalized): {rmse}");

        assert!(
            rmse < 0.0002,
            "MS Amanda score RMSE {rmse} >= 0.0002 (see msamanda_comparison.tsv)"
        );
    }
}
