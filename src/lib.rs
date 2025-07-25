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

#[cfg(test)]
mod tests {
    use std::env;

    use ndarray_stats::DeviationExt;
    use polars::prelude::*;
    use rayon::prelude::*;

    use crate::{
        configuration::Configuration,
        fast_xcorr::FastXcorr,
        utils::tests::{get_spectrum, read_test_data},
        xcorr::Xcorr,
    };

    // Test xcorr implementations agains high-res MS data
    #[test]
    fn test_xcorr() {
        let comet_df = read_test_data();

        #[allow(clippy::type_complexity)]
        let (scan_col, (peptide_col, (comet_xcorr_col, (xcorrrs_col, fast_xcorrrs_col)))): (
            Vec<i64>,
            (Vec<String>, (Vec<f64>,(Vec<f64>, Vec<f64>))),
        ) = (0..comet_df.height())
            .into_par_iter()
            .map(|idx| {
                let scan = comet_df["scan"].i64().unwrap().get(idx).unwrap();
                let comet_xcorr = comet_df["xcorr"].f64().unwrap().get(idx).unwrap();
                let proforma_peptide = comet_df["proforma_peptide"]
                    .str()
                    .unwrap()
                    .get(idx)
                    .unwrap();
                let charge = comet_df["charge"].i64().unwrap().get(idx).unwrap() as usize;

                let (mz_array, intensity_array) = get_spectrum(scan.to_string().as_str());

                let config = Configuration {
                    use_flanking_peaks: true,
                    max_fragment_charge: 5,
                    ..Configuration::default()
                };

                // xcorr implementation
                let xcorr = Xcorr::new(&config, (&mz_array, &intensity_array), charge).unwrap();
                let scoring = xcorr.xcorr_peptide(proforma_peptide).unwrap();

                // fast xcorr implementation
                let fast_xcorr =
                    FastXcorr::new(&config, (&mz_array, &intensity_array), charge).unwrap();

                let fast_scoring = fast_xcorr.xcorr_peptide(proforma_peptide).unwrap();

                (
                    scan,
                    (
                        proforma_peptide.to_string(),
                        (comet_xcorr, (scoring.score, fast_scoring.score)),
                    ),
                )
            })
            .unzip();

        let mut xcorrrs_df = DataFrame::new(vec![
            Column::new("scan".into(), scan_col),
            Column::new("modified_peptide".into(), peptide_col),
            Column::new("comet_xcorr".into(), comet_xcorr_col),
            Column::new("xcorrrs".into(), xcorrrs_col),
            Column::new("fast_xcorrrs".into(), fast_xcorrrs_col),
        ])
        .unwrap();

        CsvWriter::new(std::fs::File::create("comparison.tsv").unwrap())
            .with_separator(b'\t')
            .finish(&mut xcorrrs_df)
            .unwrap();

        if env::var("VERBOSE").is_ok() {
            let max_comet_xcorr = xcorrrs_df["comet_xcorr"].f64().unwrap().max().unwrap();
            for col_name in ["xcorrrs", "fast_xcorrrs"] {
                // Plot the comet xcorrs vs xcorrrs

                let max_calculated_xcorr = xcorrrs_df[col_name].f64().unwrap().max().unwrap();

                let mut plot = plotly::Plot::new();
                let diagonal_trace = plotly::Scatter::new(
                    vec![0.0, max_comet_xcorr],
                    vec![0.0, max_calculated_xcorr],
                )
                .mode(plotly::common::Mode::Lines)
                .marker(plotly::common::Marker::default().color("red"))
                .hover_info(plotly::common::HoverInfo::None)
                .show_legend(false);

                let correlation_trace = plotly::Scatter::new(
                    xcorrrs_df["comet_xcorr"].f64().unwrap().to_vec(),
                    xcorrrs_df[col_name].f64().unwrap().to_vec(),
                )
                .mode(plotly::common::Mode::Markers)
                .marker(plotly::common::Marker::default().color("blue"))
                .show_legend(false);

                plot.add_trace(diagonal_trace);
                plot.add_trace(correlation_trace);

                plot.set_layout(
                    plotly::Layout::new()
                        .title(format!("Comet xcorr vs {col_name}"))
                        .x_axis(
                            plotly::layout::Axis::new()
                                .title("Comet xcorr")
                                .constrain(plotly::layout::AxisConstrain::Domain),
                        )
                        .y_axis(
                            plotly::layout::Axis::new()
                                .title(col_name)
                                .scale_anchor("x"),
                        ),
                );
                plot.write_html(format!("99-{col_name}_vs_comet_xcorr.html"));
            }
        }

        // Normalize comet xcorrs and calculates xcorrs

        let comet_xcorr_max = xcorrrs_df
            .column("comet_xcorr")
            .unwrap()
            .f64()
            .unwrap()
            .max()
            .unwrap();

        let xcorrrs_max = xcorrrs_df
            .column("xcorrrs")
            .unwrap()
            .f64()
            .unwrap()
            .max()
            .unwrap();

        let fast_xcorrrs_max = xcorrrs_df
            .column("fast_xcorrrs")
            .unwrap()
            .f64()
            .unwrap()
            .max()
            .unwrap();

        let max_score = comet_xcorr_max.max(xcorrrs_max).max(fast_xcorrrs_max);

        let scaled_comet_xcorr = xcorrrs_df
            .column("comet_xcorr")
            .unwrap()
            .f64()
            .unwrap()
            .to_ndarray()
            .unwrap()
            .mapv(|x| x / max_score);

        let scaled_xcorrrs = xcorrrs_df
            .column("xcorrrs")
            .unwrap()
            .f64()
            .unwrap()
            .to_ndarray()
            .unwrap()
            .mapv(|x| x / max_score);

        let scaled_fast_xcorrrs = xcorrrs_df
            .column("fast_xcorrrs")
            .unwrap()
            .f64()
            .unwrap()
            .to_ndarray()
            .unwrap()
            .mapv(|x| x / max_score);

        let rmse_xcorrs = scaled_comet_xcorr.mean_sq_err(&scaled_xcorrrs).unwrap();
        let rmse_fast_xcorrs = scaled_comet_xcorr
            .mean_sq_err(&scaled_fast_xcorrrs)
            .unwrap();

        println!("RMSE comet xcorr vs xcorrrs: {rmse_xcorrs}");
        println!("RMSE comet xcorr vs fast xcorrrs: {rmse_fast_xcorrs}");

        assert!(rmse_xcorrs < 0.0002, "xcorr RMSE {rmse_xcorrs} >= 0.0002");
        assert!(
            rmse_fast_xcorrs < 0.0002,
            "fast xcorr RMSE {rmse_fast_xcorrs} >= 0.0002"
        );
    }
}
