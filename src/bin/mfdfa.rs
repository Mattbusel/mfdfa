//! `mfdfa` command line tool (feature `cli`).
//!
//! Reads one numeric column from a CSV file and prints the MFDFA spectrum,
//! the Higuchi fractal dimension and the wavelet Hurst exponent.
//!
//! ```text
//! mfdfa prices.csv --column close --prices
//! mfdfa returns.csv --column 0 --order 2 --q-min -4 --q-max 4 --q-step 1
//! ```

use std::process::ExitCode;

use clap::Parser;
use mfdfa::{default_scales, dfa_hurst, higuchi_fd, mfdfa, wavelet_hurst};

#[derive(Parser, Debug)]
#[command(
    name = "mfdfa",
    version,
    about = "Multifractal DFA, Higuchi dimension and wavelet Hurst for a CSV column"
)]
struct Args {
    /// CSV file ("-" for stdin).
    path: String,
    /// Column name, or zero-based index.
    #[arg(long, default_value = "0")]
    column: String,
    /// The column holds prices: analyse log returns (Higuchi uses log prices).
    #[arg(long)]
    prices: bool,
    /// The file has no header row.
    #[arg(long)]
    no_header: bool,
    /// Detrending polynomial order.
    #[arg(long, default_value_t = 1)]
    order: usize,
    /// Smallest q.
    #[arg(long, default_value_t = -5.0, allow_hyphen_values = true)]
    q_min: f64,
    /// Largest q.
    #[arg(long, default_value_t = 5.0, allow_hyphen_values = true)]
    q_max: f64,
    /// q step.
    #[arg(long, default_value_t = 0.5)]
    q_step: f64,
    /// Higuchi k_max.
    #[arg(long, default_value_t = 10)]
    k_max: usize,
}

fn read_column(args: &Args) -> Result<Vec<f64>, String> {
    let reader: Box<dyn std::io::Read> = if args.path == "-" {
        Box::new(std::io::stdin())
    } else {
        Box::new(std::fs::File::open(&args.path).map_err(|e| format!("{}: {e}", args.path))?)
    };
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(!args.no_header)
        .flexible(true)
        .from_reader(reader);
    let idx = match args.column.parse::<usize>() {
        Ok(i) => i,
        Err(_) => {
            let headers = rdr.headers().map_err(|e| e.to_string())?;
            headers
                .iter()
                .position(|h| h.trim() == args.column)
                .ok_or_else(|| format!("column '{}' not found", args.column))?
        }
    };
    let mut out = Vec::new();
    for rec in rdr.records() {
        let rec = rec.map_err(|e| e.to_string())?;
        if let Some(v) = rec.get(idx).and_then(|f| f.trim().parse::<f64>().ok()) {
            out.push(v);
        }
    }
    Ok(out)
}

fn main() -> ExitCode {
    let args = Args::parse();
    let raw = match read_column(&args) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };
    let (increments, path): (Vec<f64>, Vec<f64>) = if args.prices {
        let inc = raw.windows(2).map(|w| (w[1] / w[0]).ln()).collect();
        (inc, raw.iter().map(|p| p.ln()).collect())
    } else {
        let mut acc = 0.0;
        let path = raw
            .iter()
            .map(|v| {
                acc += v;
                acc
            })
            .collect();
        (raw, path)
    };
    if args.q_step <= 0.0 || args.q_max <= args.q_min {
        eprintln!("error: need q_max > q_min and q_step > 0");
        return ExitCode::FAILURE;
    }
    let steps = ((args.q_max - args.q_min) / args.q_step).round() as usize;
    let q: Vec<f64> = (0..=steps)
        .map(|i| args.q_min + args.q_step * i as f64)
        .collect();
    let scales = default_scales(increments.len(), args.order);

    println!(
        "n = {} increments, order = {}, scales = {:?}",
        increments.len(),
        args.order,
        scales
    );
    match mfdfa(&increments, &q, &scales, args.order) {
        Ok(spec) => {
            println!(
                "{:>8} {:>10} {:>10} {:>10} {:>10}",
                "q", "h(q)", "tau(q)", "alpha", "f(alpha)"
            );
            for i in 0..spec.q.len() {
                println!(
                    "{:>8.2} {:>10.5} {:>10.5} {:>10.5} {:>10.5}",
                    spec.q[i], spec.h_q[i], spec.tau_q[i], spec.alpha[i], spec.f_alpha[i]
                );
            }
            println!(
                "width = {:.5}  delta_h = {:.5}",
                spec.width(),
                spec.delta_h()
            );
        }
        Err(e) => println!("mfdfa: {e}"),
    }
    match dfa_hurst(&increments, &scales, args.order) {
        Ok(h) => println!("dfa_hurst h(2) = {h:.5}"),
        Err(e) => println!("dfa_hurst: {e}"),
    }
    println!(
        "higuchi_fd(k_max={}) = {:.5}",
        args.k_max,
        higuchi_fd(&path, args.k_max)
    );
    match wavelet_hurst(&increments) {
        Ok(h) => println!("wavelet_hurst = {h:.5}"),
        Err(e) => println!("wavelet_hurst: {e}"),
    }
    ExitCode::SUCCESS
}
