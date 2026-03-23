use clap::{Parser, ValueEnum};

const ASCII_ART: &str = r#"
 /$$$$$$$$ /$$$$$$$                  /$$
| $$_____/| $$__  $$                | $$
| $$      | $$  \ $$  /$$$$$$   /$$$$$$$ /$$   /$$  /$$$$$$
| $$$$$   | $$  | $$ /$$__  $$ /$$__  $$| $$  | $$ /$$__  $$
| $$__/   | $$  | $$| $$$$$$$$| $$  | $$| $$  | $$| $$  \ $$
| $$      | $$  | $$| $$_____/| $$  | $$| $$  | $$| $$  | $$
| $$      | $$$$$$$/|  $$$$$$$|  $$$$$$$|  $$$$$$/| $$$$$$$/
|__/      |_______/  \_______/ \_______/ \______/ | $$____/
                                                  | $$
                                                  | $$
                                                  |__/
"#;

/// Command Line Interface (CLI) configuration for FDedup.
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "A fast and memory-efficient FASTX PCR deduplication tool (Supports Paired-End)",
    before_help = ASCII_ART,
    arg_required_else_help = true
)]
#[command(help_expected = true)]
pub struct Cli {
    /// Path to the input FASTX file (R1 or Single-End)
    #[arg(required = true, short = '1', long)]
    pub input: String,

    /// Path to the R2 input FASTX file (Optional, enables Paired-End mode)
    #[arg(short = '2', long)]
    pub input_r2: Option<String>,

    /// Path to the output file (R1 or Single-End)
    /// Supported formats: .fastq, .fq, .fasta, .fa, .fna (+ .gz for compression)
    #[arg(short = 'o', long, default_value = "output_R1.fastq.gz")]
    pub output: String,

    /// Path to the R2 output file (Required if --input-r2 is provided)
    /// Must have the same format as --output
    #[arg(short = 'p', long)]
    pub output_r2: Option<String>,

    /// Force overwrite of output files if they exist
    #[arg(long, short)]
    pub force: bool,

    /// Enable verbose logs
    #[arg(long, short)]
    pub verbose: bool,

    /// Calculate the duplication rate without creating output files
    #[arg(long, short = 's')]
    pub dry_run: bool,

    /// Threshold for automatic hashing size selection (ignored if --hash is set)
    #[arg(long, short = 't', default_value_t = 0.001)]
    pub threshold: f64,

    /// Manually specify the hashing size (64 or 128 bits)
    #[arg(long, short = 'H')]
    pub hash: Option<HashMode>,

    /// Set the GZIP compression level (1-9)
    #[arg(long, short = 'c', default_value_t = 6, value_parser = clap::value_parser!(u32).range(1..=9))]
    pub compression: u32,
}

/// The selected hash mode (64-bit or 128-bit)
#[derive(ValueEnum, Clone, Debug)]
pub enum HashMode {
    #[value(name = "64")]
    Bit64,
    #[value(name = "128")]
    Bit128,
}