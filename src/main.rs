use csv::ReaderBuilder;
use futures::{future, StreamExt};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use reqwest::Client;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::sync::Arc;
use taxonomy::ncbi::load;
use taxonomy::Taxonomy;
use tokio::sync::Semaphore;
use tokio::task;

const ASSEMBLY_SUMMARY_PATH: &str = "assembly_summary_refseq.txt";
const TAXDUMP_DIR: &str = "taxdump";

const TARGET_TAX_ID: &str = "821"; // Phocaeicola dorei

#[derive(Debug, serde::Deserialize)]
struct NCBIAssembly {
    taxid: String,
    ftp_path: String,
    asm_name: String,
}

type BoxedError = Box<dyn Error + Send + Sync + 'static>;

async fn download_assembly(
    client: &Client,
    assembly: &NCBIAssembly,
    pb: ProgressBar,
) -> Result<(), BoxedError> {
    pb.set_message(format!("Starting {}", assembly.asm_name));

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

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let assembly_summary_file = File::open(ASSEMBLY_SUMMARY_PATH)?;

    // this seems to use the least amount of memory out of all the formats
    println!("Loading taxonomy from {}...", TAXDUMP_DIR);
    let tax = load(TAXDUMP_DIR)?;

    let descendant_tax_ids = tax.descendants(TARGET_TAX_ID)?;

    println!(
        "Found {} descendants of {} ({})...",
        descendant_tax_ids.len(),
        tax.name(TARGET_TAX_ID)?,
        TARGET_TAX_ID
    );

    println!("Reading assembly summaries from {}", ASSEMBLY_SUMMARY_PATH,);

    // skip first line because it doesn't contain an actual header
    let mut buf_reader = BufReader::new(assembly_summary_file);
    let mut first_line = String::new();
    buf_reader.read_line(&mut first_line)?;

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(buf_reader);

    let mut assemblies: Vec<NCBIAssembly> = Vec::new();

    // todo: start downloading Assemblies immediately
    for result in reader.deserialize() {
        let assembly: NCBIAssembly = result?;
        if descendant_tax_ids.contains(&assembly.taxid.as_str()) {
            let name = tax.name(assembly.taxid.as_str())?;
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
                ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{bar:.cyan/blue}] {bytes}/{total_bytes} - {msg}")
                    .unwrap()
                    .progress_chars("#>-"),
            );
            (assembly, pb)
        })
        .collect();

    // Download assemblies in parallel
    let client = Client::new();

    let semaphore = Arc::new(Semaphore::new(12));

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
