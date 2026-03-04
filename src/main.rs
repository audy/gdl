use arrow::array::StringArray;
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use clap::{ArgGroup, Parser, ValueEnum};
use csv::ReaderBuilder;
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use reqwest::blocking::Client;
use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread;
use std::time::Duration;
use tar::Archive;
use taxonomy::ncbi;
use taxonomy::parquet as tax_parquet;
use taxonomy::{GeneralTaxonomy, Taxonomy};

const TAXDUMP_URL: &str = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz";

const PB_DOWNLOAD_TEMPLATE: &str =
    "[{elapsed:.cyan}] {msg} [{bar:.green}] {bytes:.blue}/{total_bytes:.blue}";
const PB_PROGRESS_TEMPLATE: &str =
    "[{elapsed:.cyan}] {msg} [{bar:.green}] {percent:.blue}% (ETA: {eta})";
const PB_SPINNER_TEMPLATE: &str = "[{elapsed:.cyan}] {msg}";
const PROGRESS_CHARS: &str = "█░ ";

#[derive(Parser, Debug)]
#[command(group(
        ArgGroup::new("tax_id_or_name")
        .required(true)
        .args(&["tax_id", "tax_name"])
))]
struct Args {
    /// path to extracted taxdump.tar.gz
    #[clap(long, default_value = "taxdump")]
    taxdump_path: String,

    /// do not actually download anything
    #[clap(long, default_value = "false")]
    dry_run: bool,

    /// re-fetch assembly_summary.txt and taxdump
    #[clap(long, default_value = "false")]
    no_cache: bool,

    #[clap(long, default_value = "1")]
    parallel: usize,

    /*
    OUTPUT OPTIONS
    */
    #[clap(value_enum, long, default_value_t = AssemblyFormat::Fna)]
    format: AssemblyFormat,

    /// output directory, default=pwd
    #[clap(long)]
    out_dir: Option<String>,

    /*
    FILTERING PARAMETERS
    */
    /// where to fetch assemblies from (default is RefSeq)
    #[clap(value_enum, long, default_value_t = AssemblySource::Refseq)]
    source: AssemblySource,

    /// path to assembly_summary.txt
    #[clap(long)]
    assembly_summary_path: Option<String>,

    /// tax_id to download assemblies for (includes descendants unless --no-children is enabled)
    #[clap(long)]
    tax_id: Option<String>, // should this be an int (for validation)

    /// do not include child tax IDs of --tax-id (only download assemblies that have the same tax
    /// ID as provided by --tax-id)
    #[clap(long, default_value = "false")]
    no_children: bool,

    /// tax_name to download assemblies for (includes descendants unless --no-children is enabled)
    #[clap(long)]
    tax_name: Option<String>,

    /// include assemblies that match this assembly level. By default, all assembly_levels are
    /// included
    #[clap(long)]
    assembly_level: Option<Vec<String>>,
}

#[derive(Debug, serde::Deserialize)]
struct NCBIAssembly {
    taxid: String,
    ftp_path: String,
    // asm_name: String,
    assembly_level: String,
}

#[derive(ValueEnum, Clone, Debug)]
#[clap(rename_all = "lowercase")]
enum AssemblyFormat {
    Fna,
    Faa,
    Gbff,
    Gff,
}

impl AssemblyFormat {
    fn as_str(&self) -> &'static str {
        match self {
            AssemblyFormat::Fna => "fna",
            AssemblyFormat::Faa => "faa",
            AssemblyFormat::Gbff => "gbff",
            AssemblyFormat::Gff => "gff",
        }
    }
}

#[derive(ValueEnum, Clone, Debug)]
enum AssemblySource {
    Genbank,
    Refseq,
    None,
}

impl AssemblySource {
    fn as_str(&self) -> &'static str {
        match self {
            AssemblySource::Genbank => "genbank",
            AssemblySource::Refseq => "refseq",
            _ => unreachable!(),
        }
    }

    fn url(&self) -> &'static str {
        match self {
            AssemblySource::Genbank => {
                "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
            }
            AssemblySource::Refseq => {
                "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"
            }
            _ => unreachable!(),
        }
    }
}

// here we should re-use a single client to take advantage of keep-alive connection pooling
fn download_assembly(
    client: &Client,
    assembly: &NCBIAssembly,
    format: &AssemblyFormat,
    out_path: &Path,
) -> PathBuf {
    // TODO: use a proper url parser
    let last_part = assembly.ftp_path.split('/').last().unwrap_or_else(|| {
        panic!(
            "Failed to get the filename from FTP path {}",
            assembly.ftp_path
        )
    });

    let url = format!(
        "{}/{}_genomic.{}.gz",
        assembly.ftp_path,
        last_part,
        format.as_str()
    );

    let assembly_filename = format!("{}.{}.gz", last_part, format.as_str());
    let assembly_path = out_path.join(assembly_filename);

    let mut file = File::create(&assembly_path)
        .unwrap_or_else(|_| panic!("Unable to write to {}", assembly_path.display()));

    let mut response = client
        .get(&url)
        .send()
        .unwrap_or_else(|_| panic!("Error fetching data from {}", url));

    response
        .copy_to(&mut file)
        .unwrap_or_else(|_| panic!("Unable to write to {}", assembly_path.display()));

    assembly_path
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
                1 => Ok(matches
                    .first()
                    .unwrap_or_else(|| panic!("No tax ID found for name {}", tax_name))),
                // TODO: show matched lineages and their tax IDs to help the user disambiguate
                _ => Err("Name is ambiguous"),
            }
        }
        _ => Err("Either --tax-id or --tax-name must be provided, but not both"),
    }
}

fn download_and_extract_taxdump(path: &str) {
    let client = Client::new();
    let mut response = client
        .get(TAXDUMP_URL)
        .send()
        .unwrap_or_else(|_| panic!("Unable to fetch NCBI taxonomy dump from {}", TAXDUMP_URL));

    let content_length = response.content_length().unwrap_or(0);

    let pb = ProgressBar::new(content_length);
    pb.set_style(
        ProgressStyle::with_template(PB_DOWNLOAD_TEMPLATE)
            .unwrap()
            .progress_chars(PROGRESS_CHARS),
    );
    pb.set_message("taxdump.tar.gz");

    let file = File::create("taxdump.tar.gz").expect("Unable to read taxdump.tar.gz");
    let mut wrapped_file = pb.wrap_write(file);

    let _ = response.copy_to(&mut wrapped_file);

    pb.set_message("Extracting taxonomy");
    let tar_gz = File::open("taxdump.tar.gz").expect("Unable to open taxdump.tar.gz");
    let decompressed = GzDecoder::new(tar_gz);
    let mut archive = Archive::new(decompressed);

    std::fs::create_dir_all(path)
        .unwrap_or_else(|_| panic!("Unable to create taxdump output dir: {}", path));
    archive
        .unpack(path)
        .expect("Unable to extract taxdump.tar.gz");

    fs::remove_file("taxdump.tar.gz").expect("Unable to remove taxdump.tar.gz");

    pb.set_message("Converting taxonomy to parquet");
    let tax = ncbi::load(path).unwrap_or_else(|_| panic!("Unable to load taxdump from {}", path));
    let parquet_path = format!("{}.parquet", path);
    tax_parquet::save::<&str, _, _>(&tax, &parquet_path)
        .unwrap_or_else(|_| panic!("Unable to save taxonomy to {}", parquet_path));
    fs::remove_dir_all(path)
        .unwrap_or_else(|_| panic!("Unable to remove taxdump directory {}", path));

    pb.finish();
}

fn download_assembly_summary(assembly_source: &AssemblySource, out_path: &str) {
    let client = Client::new();

    let assembly_summary_url = assembly_source.url();

    let mut response = client.get(assembly_summary_url).send().unwrap_or_else(|_| {
        panic!(
            "Unable to fetch assembly summary from {}",
            assembly_summary_url
        )
    });

    let content_length = response.content_length().unwrap_or(0);

    let pb = ProgressBar::new(content_length);
    pb.set_style(
        ProgressStyle::with_template(PB_DOWNLOAD_TEMPLATE)
            .unwrap()
            .progress_chars(PROGRESS_CHARS),
    );

    pb.set_message(out_path.to_string());

    let file = File::create(out_path)
        .unwrap_or_else(|_| panic!("Unable to open assembly summary {}", out_path));
    let mut wrapped_file = pb.wrap_write(file);

    let _ = response.copy_to(&mut wrapped_file);

    pb.finish();
}

fn assembly_summary_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("taxid", DataType::Utf8, false),
        Field::new("ftp_path", DataType::Utf8, false),
        Field::new("assembly_level", DataType::Utf8, false),
    ]))
}

fn save_assembly_summary_to_parquet(txt_path: &str, parquet_path: &str) {
    let assembly_summary_file = File::open(txt_path)
        .unwrap_or_else(|_| panic!("Unable to open assembly summary {}", txt_path));

    let mut buf_reader = BufReader::new(assembly_summary_file);
    let mut first_line = String::new();
    buf_reader
        .read_line(&mut first_line)
        .expect("Unable to skip first line of assembly summary");

    let pb = ProgressBar::new(
        buf_reader
            .get_ref()
            .metadata()
            .expect("Unable to get file size")
            .len(),
    );
    pb.set_style(
        ProgressStyle::with_template(PB_PROGRESS_TEMPLATE)
            .unwrap()
            .progress_chars(PROGRESS_CHARS),
    );
    pb.set_message(format!("Converting {} to parquet", txt_path));

    let wrapped_reader = pb.wrap_read(buf_reader);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(wrapped_reader);

    let mut taxids: Vec<String> = Vec::new();
    let mut ftp_paths: Vec<String> = Vec::new();
    let mut assembly_levels: Vec<String> = Vec::new();

    for result in reader.deserialize() {
        let assembly: NCBIAssembly = result.expect("Unable to parse assembly summary line");
        taxids.push(assembly.taxid);
        ftp_paths.push(assembly.ftp_path);
        assembly_levels.push(assembly.assembly_level);
    }

    let schema = assembly_summary_schema();
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(taxids)),
            Arc::new(StringArray::from(ftp_paths)),
            Arc::new(StringArray::from(assembly_levels)),
        ],
    )
    .expect("Unable to create assembly summary record batch");

    let file = File::create(parquet_path)
        .unwrap_or_else(|_| panic!("Unable to create parquet file {}", parquet_path));
    let mut writer =
        ArrowWriter::try_new(file, schema, None).expect("Unable to create ArrowWriter");
    writer
        .write(&batch)
        .expect("Unable to write assembly summary parquet");
    writer
        .close()
        .expect("Unable to close assembly summary parquet writer");

    pb.finish_with_message(format!("Saved assembly summary to {}", parquet_path));
}

fn filter_assemblies_from_parquet(
    parquet_path: &str,
    filter_assembly_levels: Option<Vec<String>>,
    filter_tax_ids: HashSet<&str>,
) -> Vec<NCBIAssembly> {
    let file = File::open(parquet_path)
        .unwrap_or_else(|_| panic!("Unable to open {}", parquet_path));
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .unwrap_or_else(|_| panic!("Unable to read parquet file {}", parquet_path));
    let total_rows = builder.metadata().file_metadata().num_rows();
    let reader = builder
        .build()
        .expect("Unable to build parquet record batch reader");

    let pb = ProgressBar::new(total_rows as u64);
    pb.set_style(
        ProgressStyle::with_template(PB_PROGRESS_TEMPLATE)
            .unwrap()
            .progress_chars(PROGRESS_CHARS),
    );
    pb.set_message(format!("Filtering {}", parquet_path));

    let mut assemblies: Vec<NCBIAssembly> = Vec::new();

    for result in reader {
        let batch = result.expect("Unable to read parquet batch");
        let taxid_col = batch
            .column_by_name("taxid")
            .expect("Missing 'taxid' column")
            .as_any()
            .downcast_ref::<StringArray>()
            .expect("'taxid' must be Utf8");
        let ftp_path_col = batch
            .column_by_name("ftp_path")
            .expect("Missing 'ftp_path' column")
            .as_any()
            .downcast_ref::<StringArray>()
            .expect("'ftp_path' must be Utf8");
        let assembly_level_col = batch
            .column_by_name("assembly_level")
            .expect("Missing 'assembly_level' column")
            .as_any()
            .downcast_ref::<StringArray>()
            .expect("'assembly_level' must be Utf8");

        for i in 0..batch.num_rows() {
            let taxid = taxid_col.value(i);
            let assembly_level = assembly_level_col.value(i);
            if filter_tax_ids.contains(taxid)
                && (filter_assembly_levels.is_none()
                    || filter_assembly_levels
                        .as_ref()
                        .unwrap()
                        .contains(&assembly_level.to_string()))
            {
                assemblies.push(NCBIAssembly {
                    taxid: taxid.to_string(),
                    ftp_path: ftp_path_col.value(i).to_string(),
                    assembly_level: assembly_level.to_string(),
                });
            }
        }

        pb.inc(batch.num_rows() as u64);
    }

    pb.finish_with_message(format!("Kept {} assemblies", assemblies.len()));

    assemblies
}

fn load_taxonomy(taxdump_path: &str) -> GeneralTaxonomy {
    let parquet_path = format!("{}.parquet", taxdump_path);
    if Path::new(&parquet_path).exists() {
        tax_parquet::load(&parquet_path)
            .unwrap_or_else(|_| panic!("Unable to load taxonomy from {}", parquet_path))
    } else {
        ncbi::load(taxdump_path)
            .unwrap_or_else(|_| panic!("Unable to load taxdump from {}", taxdump_path))
    }
}

fn filter_assemblies(
    assembly_summary_path: &String,
    filter_assembly_levels: Option<Vec<String>>,
    filter_tax_ids: HashSet<&str>,
) -> Vec<NCBIAssembly> {
    if assembly_summary_path.ends_with(".parquet") {
        filter_assemblies_from_parquet(assembly_summary_path, filter_assembly_levels, filter_tax_ids)
    } else {
        filter_assemblies_from_txt(assembly_summary_path, filter_assembly_levels, filter_tax_ids)
    }
}

fn filter_assemblies_from_txt(
    assembly_summary_path: &String,
    // TODO: combine multiple with AND/OR?
    filter_assembly_levels: Option<Vec<String>>,
    filter_tax_ids: HashSet<&str>,
) -> Vec<NCBIAssembly> {
    // filter assembly summaries
    let assembly_summary_file = File::open(&assembly_summary_path).unwrap_or_else(|_| {
        panic!(
            "Unable to open assembly summary path {}",
            assembly_summary_path
        )
    });

    // skip first line because it doesn't contain an actual header
    let mut buf_reader = BufReader::new(assembly_summary_file);
    let mut first_line = String::new();

    buf_reader
        .read_line(&mut first_line)
        .expect("Unable to parse assembly summaries");

    let pb = ProgressBar::new(
        buf_reader
            .get_ref()
            .metadata()
            .expect("Unable to get file size")
            .len(),
    );
    pb.set_style(
        ProgressStyle::with_template(PB_PROGRESS_TEMPLATE)
            .unwrap()
            .progress_chars(PROGRESS_CHARS),
    );
    pb.set_message(format!("Filtering {}", assembly_summary_path));

    let wrapped_reader = pb.wrap_read(buf_reader);

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(wrapped_reader);

    let mut assemblies: Vec<NCBIAssembly> = Vec::new();

    for result in reader.deserialize() {
        let assembly: NCBIAssembly = result.expect("Unable to parse assembly summary line");

        if filter_tax_ids.contains(&assembly.taxid.as_str())
            && (filter_assembly_levels.is_none()
                || (filter_assembly_levels
                    .as_ref()
                    .expect("Unable to parse assembly level")
                    .contains(&assembly.assembly_level)))
        {
            assemblies.push(assembly);
        }
    }

    pb.finish_with_message(format!("Kept {} assemblies", assemblies.len()));

    assemblies
}

fn main() {
    let args = Args::parse();

    // either use the provided assembly summary file or fetch it from source. if fetching from
    // source and it already exists; just use the cached parquet unless --no-cache is enabled.
    let assembly_summary_path = match (args.assembly_summary_path, &args.source) {
        (None, assembly_source) => {
            let txt_path = format!("assembly_summary_{}.txt", assembly_source.as_str());
            let parquet_path = format!("assembly_summary_{}.parquet", assembly_source.as_str());
            if args.no_cache || !Path::new(&parquet_path).exists() {
                download_assembly_summary(assembly_source, &txt_path);
                save_assembly_summary_to_parquet(&txt_path, &parquet_path);
                fs::remove_file(&txt_path)
                    .unwrap_or_else(|_| panic!("Unable to remove {}", txt_path));
            };
            parquet_path
        }
        (Some(assembly_summary_path), AssemblySource::None) => assembly_summary_path,
        _ => {
            panic!("--source and --assembly-summary-path are mutually exclusive")
        }
    };

    // download taxonomy
    let taxdump_parquet_path = format!("{}.parquet", &args.taxdump_path);
    if args.no_cache
        || (!Path::new(&args.taxdump_path).exists() && !Path::new(&taxdump_parquet_path).exists())
    {
        download_and_extract_taxdump(&args.taxdump_path);
    }

    let pb = ProgressBar::new(0);
    pb.set_style(ProgressStyle::with_template(PB_SPINNER_TEMPLATE).unwrap());
    pb.set_message(format!("Loading taxonomy from {}", &args.taxdump_path));

    // Spawn a separate thread to tick the spinner
    let pb_clone = pb.clone();
    thread::spawn(move || {
        while !pb_clone.is_finished() {
            pb_clone.tick();
            thread::sleep(Duration::from_millis(100));
        }
    });

    let tax = load_taxonomy(&args.taxdump_path);

    let tax_id: &str = get_tax_id(args.tax_id.as_deref(), args.tax_name.as_deref(), &tax)
        .expect("Unable to find a tax ID");

    pb.finish_with_message(format!("Loaded {} taxa", tax.names.len()));

    let descendant_tax_ids: HashSet<&str> = if args.no_children {
        [tax_id].into()
    } else {
        tax.descendants(tax_id)
            .unwrap_or_else(|_| {
                panic!("Unable to find taxonomic descendants for tax ID {}", tax_id)
            })
            .into_iter()
            .chain([tax_id])
            .collect()
    };

    let assemblies = filter_assemblies(
        &assembly_summary_path,
        args.assembly_level,
        descendant_tax_ids,
    );

    let n_assemblies = assemblies.len();

    // setup threadpool using --parallel
    ThreadPoolBuilder::new()
        .num_threads(args.parallel)
        .build_global()
        .expect("Unable to build thread pool");

    let out_dir = args.out_dir.unwrap_or(".".to_string());
    let out_path = Path::new(&out_dir);

    if !out_path.exists() {
        fs::create_dir_all(out_path).expect("Unable to create path");
    }

    if !args.dry_run {
        // Download assemblies in parallel
        let client = Client::new();

        let pb = ProgressBar::new(n_assemblies as u64);
        pb.set_style(
            ProgressStyle::with_template(PB_PROGRESS_TEMPLATE)
                .unwrap()
                .progress_chars(PROGRESS_CHARS),
        );
        pb.set_message(format!(
            "Downloading {} assemblies within the {} `{}` (tax_id={}) in {} format\n",
            assemblies.len(),
            tax.rank(tax_id).unwrap(),
            tax.name(tax_id).unwrap(),
            tax_id,
            &args.format.as_str()
        ));
        let _tasks: Vec<_> = assemblies
            .par_iter()
            .map(|assembly| {
                let client = client.clone();
                pb.inc(1);
                let _ = download_assembly(&client, assembly, &args.format, out_path);
            })
            .collect();

        pb.finish_with_message(format!(
            "Saved {} assemblies to {}",
            assemblies.len(),
            out_dir
        ));
    }

    println!("Thank you for flying gdl!");
}

#[cfg(test)]
mod tests {
    use super::*;
    use httpmock::Method::GET;
    use httpmock::MockServer;
    use tempfile::tempdir;

    #[test]
    fn test_download_assembly() {
        // Start a mock HTTP server
        let server = MockServer::start();
        let file_content = b"test genome data";
        let ftp_path = format!("{}/test_asm", server.url(""));
        let mock = server.mock(|when, then| {
            when.method(GET).path("/test_asm/test_asm_genomic.fna.gz");
            then.status(200)
                .header("Content-Type", "application/octet-stream")
                .body(file_content);
        });

        let assembly = NCBIAssembly {
            taxid: "123".to_string(),
            ftp_path: ftp_path.clone(),
            assembly_level: "Complete Genome".to_string(),
        };
        let format = AssemblyFormat::Fna;
        let client = Client::new();
        let tmp_dir = tempdir().unwrap();
        let out_path = tmp_dir.path();

        let result_path = download_assembly(&client, &assembly, &format, out_path);
        let result_data = std::fs::read(&result_path).unwrap();
        assert_eq!(result_data, file_content);
        mock.assert();
    }
}
