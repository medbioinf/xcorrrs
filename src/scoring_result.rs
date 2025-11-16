use std::fmt::Display;

pub struct ScoringResult {
    pub score: f64,
    pub ions_total: usize,
    pub ions_match_tolerance_based: usize,
    pub ions_matched_sp_score_based: usize,
    pub min_theoretical_mass: f64,
    pub max_theoretical_mass: f64,
}

impl ScoringResult {
    pub fn round_score(&self, decimals: u16) -> f64 {
        let factor = 10f64.powi(decimals as i32);
        (self.score * factor).round() / factor
    }
}

impl Display for ScoringResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Score: {:.2}, Ions Total: {}, Ions Matched: {} / {} (tolerance / sp_score), Min Mass: {:.2}, Max Mass: {:.2}",
            self.score,
            self.ions_total,
            self.ions_match_tolerance_based,
            self.ions_matched_sp_score_based,
            self.min_theoretical_mass,
            self.max_theoretical_mass
        )
    }
}
