use clap::{ArgGroup, Parser};
use csv::ReaderBuilder;
use futures::{future, StreamExt};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use reqwest::Client;
use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::sync::Arc;
use taxonomy::ncbi::load;
use taxonomy::{GeneralTaxonomy, Taxonomy};
use tokio::sync::Semaphore;
use tokio::task;

#[derive(Parser, Debug)]
#[command(group(
        ArgGroup::new("tax_id_or_name")
        .required(true)
        .args(&["tax_id", "tax_name"])
))]
struct Args {
    /// path to assembly_summary.txt
    #[clap(short, long, default_value = "assembly_summary_refseq.txt")]
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
    asm_name: String,
    assembly_level: String,
}

type BoxedError = Box<dyn Error + Send + Sync + 'static>;

async fn download_assembly(
    client: &Client,
    assembly: &NCBIAssembly,
    pb: ProgressBar,
) -> Result<(), BoxedError> {
    pb.set_message(format!("Starting {}", assembly.asm_name));

    // TODO: use a proper url parser
    let last_part = assembly
        .ftp_path
        .split('/')
        .last()
        .expect("Failed to get the filename");
    let url = format!("{}/{}_genomic.fna.gz", assembly.ftp_path, last_part);

    let response = client.get(url).send().await?;
    let total_size = response.content_length().unwrap_or(0);

    pb.set_length(total_size);

    let mut file = File::create(format!("{}.fna.gz", last_part))?;
    let mut stream = response.bytes_stream();

    while let Some(chunk) = stream.next().await {
        let chunk = chunk?;
        pb.inc(chunk.len() as u64);
        file.write_all(&chunk)?;
    }

    pb.finish_and_clear();
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

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // this seems to use the least amount of memory out of all the formats
    println!("Loading taxonomy from {}...", args.taxdump_path);
    let tax = load(&args.taxdump_path)?;

    let tax_id: &str = get_tax_id(args.tax_id.as_deref(), args.tax_name.as_deref(), &tax)
        .expect("Unable to find a tax ID");

    println!("Found tax ID {}", tax_id);

    let assembly_summary_file = File::open(args.assembly_summary_path.clone())?;

    let descendant_tax_ids: HashSet<&str> = tax.descendants(tax_id)?.into_iter().collect();

    println!(
        "Found {} descendants of {} ({})...",
        descendant_tax_ids.len(),
        tax.name(tax_id)?,
        tax_id
    );

    println!(
        "Reading assembly summaries from {}",
        args.assembly_summary_path
    );

    // TODO: progress bar for ^^^

    // skip first line because it doesn't contain an actual header
    let mut buf_reader = BufReader::new(assembly_summary_file);
    let mut first_line = String::new();
    // do we really have to read it _into_ something?
    buf_reader.read_line(&mut first_line)?;

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(buf_reader);

    let mut assemblies: Vec<NCBIAssembly> = Vec::new();

    // todo: start downloading Assemblies immediately (except we do't know the total?)
    for result in reader.deserialize() {
        let assembly: NCBIAssembly = result?;
        if descendant_tax_ids.contains(&assembly.taxid.as_str())
            && (args.assembly_level.is_none()
                || (args
                    .assembly_level
                    .as_ref()
                    .expect("What")
                    .contains(&assembly.assembly_level)))
        {
            assemblies.push(assembly);
        }
    }

    println!("Downloading {} assemblies", assemblies.len());

    // TODO: one master progress bar with ETA?

    let multi_process = MultiProgress::new();
    let progress_bars: Vec<_> = assemblies
        .into_iter()
        .map(|assembly| {
            let pb = multi_process.add(ProgressBar::new(0));
            pb.set_style(
                // todo add the assembly # and total assemblies here
                ProgressStyle::with_template("{spinner:.green} [{len}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} - {msg}")
                    .unwrap()
                    .progress_chars("#>-"),
            );
            (assembly, pb)
        })
        .collect();

    // Download assemblies in parallel
    let client = Client::new();

    let semaphore = Arc::new(Semaphore::new(args.parallel as usize));

    let tasks: Vec<_> = progress_bars
        .into_iter()
        .map(|(assembly, pb)| {
            let client = client.clone();
            let semaphore = semaphore.clone();
            task::spawn(async move {
                let _permit = semaphore.acquire().await.unwrap();
                download_assembly(&client, &assembly, pb).await
            })
        })
        .collect();

    future::join_all(tasks).await;

    Ok(())
}
