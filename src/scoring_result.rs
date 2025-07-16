use std::fmt::Display;

pub struct ScoringResult {
    pub score: f64,
    pub ions_total: usize,
    pub ions_matched: usize,
    pub min_theoretical_mass: f64,
    pub max_theoretical_mass: f64,
}

impl Display for ScoringResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Score: {:.2}, Ions Total: {}, Ions Matched: {}, Min Mass: {:.2}, Max Mass: {:.2}",
            self.score,
            self.ions_total,
            self.ions_matched,
            self.min_theoretical_mass,
            self.max_theoretical_mass
        )
    }
}
