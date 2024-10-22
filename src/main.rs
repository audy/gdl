use clap::{ArgGroup, Parser};
use csv::ReaderBuilder;
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use reqwest::blocking::Client;
use std::collections::HashSet;
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::thread;
use std::time::Duration;
use tar::Archive;
use taxonomy::ncbi::load;
use taxonomy::{GeneralTaxonomy, Taxonomy};

const TAXDUMP_URL: &str = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz";
const ASSEMBLY_SUMMARY_URL: &str =
    "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt";

const PB_DOWNLOAD_TEMPLATE: &str =
    "{msg} [{elapsed_precise}] [{bar:.white/green}] {bytes}/{total_bytes}";
const PB_PROGRESS_TEMPLATE: &str = "{msg} [{elapsed_precise}] [{bar:.white/green}] {pos}/{len}";
const PB_SPINNER_TEMPLATE: &str = "{spinner:.green} {msg}";

#[derive(Parser, Debug)]
#[command(group(
        ArgGroup::new("tax_id_or_name")
        .required(true)
        .args(&["tax_id", "tax_name"])
))]
struct Args {
    /// path to assembly_summary.txt
    #[clap(short, long, default_value = "assembly_summary_genbank.txt")]
    assembly_summary_path: String,

    /// path to extracted taxdump.tar.gz
    #[clap(short, long, default_value = "taxdump")]
    taxdump_path: String,

    /// number of simultaneous downloads
    #[clap(short, long, default_value_t = 4, value_parser = clap::value_parser!(u32))]
    parallel: u32,

    /*
    FILTERING PARAMETERS
    */
    /// tax_id to download assemblies for (includes descendants)
    #[clap(short, long)]
    tax_id: Option<String>, // should this be an int (for validation)
    /// tax_name to download assemblies for (includes descendants)
    #[clap(short, long)]
    tax_name: Option<String>,

    /// include assemblies that match this assembly level. can be used multiple times
    /// by default, all assembly_levels are included
    #[clap(short, long)]
    assembly_level: Option<Vec<String>>,
}

#[derive(Debug, serde::Deserialize)]
struct NCBIAssembly {
    taxid: String,
    ftp_path: String,
    // asm_name: String,
    assembly_level: String,
}

type BoxedError = Box<dyn Error + Send + Sync + 'static>;

// here we should re-use a single client to take advantage of keep-alive connection pooling
fn download_assembly(client: &Client, assembly: &NCBIAssembly) -> Result<(), BoxedError> {
    // TODO: use a proper url parser
    let last_part = assembly
        .ftp_path
        .split('/')
        .last()
        .expect("Failed to get the filename");
    let url = format!("{}/{}_genomic.fna.gz", assembly.ftp_path, last_part);

    let mut file = File::create(format!("{}.fna.gz", last_part))?;
    let mut response = client.get(url).send()?;
    response.copy_to(&mut file)?;

    Ok(())
}

fn get_tax_id<'a>(
    tax_id: Option<&'a str>,
    tax_name: Option<&'a str>,
    tax: &'a GeneralTaxonomy,
) -> Result<&'a str, &'a str> {
    // TODO: make sure tax ID exists
    match (tax_id, tax_name) {
        (Some(tax_id), None) => Ok(tax_id),
        (None, Some(tax_name)) => {
            let matches = tax.find_all_by_name(tax_name);
            match matches.len() {
                0 => Err("No matches found"),
                1 => Ok(matches.first().expect("No tax ID?")),
                // TODO: show matched lineages and their tax IDs to help the user disambiguate
                _ => Err("Ambiguous Name!"),
            }
        }
        _ => Err("Either --tax-id or --tax-name must be provided, but not both"),
    }
}

fn download_and_extract_taxdump(path: &str) -> Result<(), BoxedError> {
    if Path::new(path).exists() {
        return Ok(());
    }

    // TODO recycle client for Keep Alive
    let client = Client::new();

    let pb = ProgressBar::new(0);
    pb.set_message("Fetching taxonomy");

    let mut response = client.get(TAXDUMP_URL).send()?;
    let mut file = File::create("taxdump.tar.gz")?;

    let _ = response.copy_to(&mut file);

    pb.set_message("Extracting taxonomy");
    let tar_gz = File::open("taxdump.tar.gz")?;
    let decompressed = GzDecoder::new(tar_gz);
    let mut archive = Archive::new(decompressed);
    std::fs::create_dir_all(path)?;
    archive.unpack(path)?;
    fs::remove_file("taxdump.tar.gz")?;

    pb.finish_and_clear();

    Ok(())
}

fn download_assembly_summary(path: &str) -> Result<(), BoxedError> {
    if Path::new(path).exists() {
        return Ok(());
    }

    // TODO: re-use existing Client
    let client = Client::new();
    let mut response = client.get(ASSEMBLY_SUMMARY_URL).send()?;

    let mut file = File::create(path)?;
    let _ = response.copy_to(&mut file);

    let pb = ProgressBar::new(0);
    pb.set_style(
        ProgressStyle::with_template(PB_DOWNLOAD_TEMPLATE)
            .unwrap()
            .progress_chars("#>-"),
    );

    pb.set_message("Fetching assembly summary");

    pb.finish_and_clear();
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // download taxonomy and assembly summary
    let _ = download_and_extract_taxdump(&args.taxdump_path);
    let _ = download_assembly_summary(&args.assembly_summary_path);

    let pb = ProgressBar::new(0);

    pb.set_style(ProgressStyle::with_template(PB_SPINNER_TEMPLATE).unwrap());
    pb.set_message("Loading taxonomy");

    // Spawn a separate thread to tick the spinner
    let pb_clone = pb.clone();
    thread::spawn(move || {
        while !pb_clone.is_finished() {
            pb_clone.tick();
            thread::sleep(Duration::from_millis(100));
        }
    });

    let tax = load(&args.taxdump_path)?;
    pb.set_message(format!("Loaded {} taxa", Taxonomy::<&str>::len(&tax)));

    let tax_id: &str = get_tax_id(args.tax_id.as_deref(), args.tax_name.as_deref(), &tax)
        .expect("Unable to find a tax ID");

    let descendant_tax_ids: HashSet<&str> = tax.descendants(tax_id)?.into_iter().collect();

    pb.finish_with_message(format!(
        "Found {} descendants of {} ({})",
        descendant_tax_ids.len(),
        tax.name(tax_id)?,
        tax_id
    ));

    // filter assembly summaries
    let assembly_summary_file = File::open(args.assembly_summary_path.clone())?;

    // TODO: progress bar for ^^^

    // skip first line because it doesn't contain an actual header
    let mut buf_reader = BufReader::new(assembly_summary_file);
    let mut first_line = String::new();
    // do we really have to read it _into_ something?
    buf_reader.read_line(&mut first_line)?;

    let pb = ProgressBar::new(buf_reader.get_ref().metadata()?.len());

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(buf_reader);

    let mut assemblies: Vec<NCBIAssembly> = Vec::new();

    pb.set_style(ProgressStyle::with_template(PB_SPINNER_TEMPLATE).unwrap());

    for result in reader.deserialize() {
        let assembly: NCBIAssembly = result?;

        if descendant_tax_ids.contains(&assembly.taxid.as_str())
            && (args.assembly_level.is_none()
                || (args
                    .assembly_level
                    .as_ref()
                    .expect("Unable to parse assembly level")
                    .contains(&assembly.assembly_level)))
        {
            assemblies.push(assembly);
            pb.set_message(format!("found {} assemblies", assemblies.len()));
        }
    }

    let n_assemblies = assemblies.len();

    pb.finish_with_message(format!("Found {n_assemblies} assemblies"));

    // Download assemblies in parallel
    let client = Client::new();

    let pb = ProgressBar::new(n_assemblies as u64);
    pb.set_style(ProgressStyle::with_template(PB_PROGRESS_TEMPLATE).unwrap());

    let _tasks: Vec<_> = assemblies
        .into_iter()
        .map(|assembly| {
            let client = client.clone();
            pb.inc(1);
            let _ = download_assembly(&client, &assembly);
        })
        .collect();

    pb.finish_with_message("Assemblies downloaded");

    Ok(())
}
