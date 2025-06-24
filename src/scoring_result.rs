pub struct ScoringResult {
    pub score: f64,
    pub ions_total: usize,
    pub ions_matched: usize,
    pub min_theoretical_mass: f64,
    pub max_theoretical_mass: f64,
}
