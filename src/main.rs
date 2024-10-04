use csv::ReaderBuilder;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use taxonomy::ncbi::load;
use taxonomy::Taxonomy;

const ASSEMBLY_SUMMARY_PATH: &str = "assembly_summary_refseq.txt";
const TAXDUMP_DIR: &str = "taxdump";

const TARGET_TAX_ID: &str = "357276"; // Phocaeicola dorei

#[derive(Debug, serde::Deserialize)]
struct NCBIAssembly {
    taxid: String,
    ftp_path: String,
}

fn main() -> Result<(), Box<dyn Error>> {
    let assembly_summary_file = File::open(ASSEMBLY_SUMMARY_PATH)?;

    let tax = load(TAXDUMP_DIR)?;

    let descendant_tax_ids = tax.descendants(TARGET_TAX_ID)?;

    println!(
        "Found {} descendants of tax ID {} ({})",
        descendant_tax_ids.len(),
        TARGET_TAX_ID,
        tax.name(TARGET_TAX_ID)?,
    );

    // skip first line because it doesn't contain an actual header
    let mut buf_reader = BufReader::new(assembly_summary_file);
    let mut first_line = String::new();
    buf_reader.read_line(&mut first_line)?;

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(buf_reader);

    for result in reader.deserialize() {
        let assembly: NCBIAssembly = result?;
        if descendant_tax_ids.contains(&assembly.taxid.as_str()) {
            let name = tax.name(assembly.taxid.as_str())?;
            println!("taxid: {} -> {}", assembly.taxid, name);
            println!("  ftp_path: {}", assembly.ftp_path);
        }
    }

    Ok(())
}
