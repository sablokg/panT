use crate::structfile::SamAlign;
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::path::Path;
use tokio::fs::File;
use tokio::io::{self, AsyncBufReadExt, AsyncWriteExt, BufReader, BufWriter};

/*
Gaurav Sablok
gsablok@proton.me
*/

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct SAMFile {
    pub alignedsame: String,
    pub referencefasta: String,
}

#[derive(Debug)]
struct SamLineError(String);

impl fmt::Display for SamLineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}
impl Error for SamLineError {}

impl SAMFile {
    pub async fn parsesam(&self) -> Result<String, Box<dyn Error>> {
        let referenceunpack = read_fasta(&self.referencefasta).await?;

        let alignmentfile = File::open(&self.alignedsame).await?;
        let mut alignmentread = BufReader::with_capacity(1 << 20, alignmentfile).lines();

        // Read every alignment line up front so the (CPU-bound) parsing
        // and sequence-slicing work below can be spread across all
        // available cores with rayon, instead of doing one line at a
        // time on a single thread interleaved with file I/O.
        let mut rawlines: Vec<String> = Vec::new();
        while let Some(line) = alignmentread.next_line().await? {
            if line.is_empty() || line.starts_with('@') {
                continue;
            }
            rawlines.push(line);
        }

        // HYPER-THREADING: turn each SAM line into its formatted
        // reference/query fasta output in parallel. `referenceunpack` is
        // only read (never mutated) from every thread, which is safe to
        // share this way. `par_iter` preserves line order, so the output
        // files come out in the same order as the input SAM file.
        let results: Vec<Result<(Option<String>, String), SamLineError>> = rawlines
            .par_iter()
            .map(|line| -> Result<(Option<String>, String), SamLineError> {
                let linevec: Vec<&str> = line.split('\t').collect();

                if linevec.len() < 10 {
                    return Err(SamLineError(format!(
                        "malformed SAM line, expected at least 10 fields: {line}"
                    )));
                }

                let pos: usize = linevec[3]
                    .parse()
                    .map_err(|e| SamLineError(format!("invalid POS field in line {line:?}: {e}")))?;
                let refname = linevec[2];
                let queryname = linevec[0];
                let queryseq = linevec[9];

                let refstart = pos.saturating_sub(1);
                let refend = refstart + queryseq.len();

                let alignment = SamAlign {
                    refname: refname.to_string(),
                    queryname: queryname.to_string(),
                    querystart: pos,
                    queryend: pos + queryseq.len(),
                    refstart,
                    refend,
                    queryseq: queryseq.to_string(),
                };

                let refrecord = referenceunpack.get(&alignment.refname).and_then(|record| {
                    record
                        .sequence
                        .get(alignment.refstart..alignment.refend)
                        .map(|seq| format!(">{}\n{}\n", alignment.refname, seq))
                });

                let queryrecord =
                    format!(">{}\n{}\n", alignment.queryname, alignment.queryseq);

                Ok((refrecord, queryrecord))
            })
            .collect();

        let refwrite_file = File::create("reference.fasta").await?;
        let querywrite_file = File::create("query.fasta").await?;
        let mut refwrite = BufWriter::with_capacity(1 << 20, refwrite_file);
        let mut querywrite = BufWriter::with_capacity(1 << 20, querywrite_file);

        for result in results {
            let (refrecord, queryrecord) = result.map_err(|e| Box::new(e) as Box<dyn Error>)?;
            if let Some(seq) = refrecord {
                refwrite.write_all(seq.as_bytes()).await?;
            }
            querywrite.write_all(queryrecord.as_bytes()).await?;
        }

        refwrite.flush().await?;
        querywrite.flush().await?;

        Ok("sam alignment has been parsed".to_string())
    }
}

#[derive(Debug)]
pub struct FastaRecord {
    pub id: String,
    pub sequence: String,
}

pub async fn read_fasta<P: AsRef<Path>>(path: P) -> io::Result<HashMap<String, FastaRecord>> {
    let file = File::open(path).await?;
    let mut reader = BufReader::with_capacity(1 << 20, file).lines();

    let mut records = HashMap::new();
    let mut current_id = String::new();
    let mut current_sequence = String::new();

    while let Some(line) = reader.next_line().await? {
        if let Some(id) = line.strip_prefix('>') {
            if !current_id.is_empty() {
                records.insert(
                    current_id.clone(),
                    FastaRecord {
                        id: std::mem::take(&mut current_id),
                        sequence: std::mem::take(&mut current_sequence),
                    },
                );
            }
            current_id = id.to_string();
        } else {
            current_sequence.push_str(&line);
        }
    }
    if !current_id.is_empty() {
        records.insert(
            current_id.clone(),
            FastaRecord {
                id: current_id,
                sequence: current_sequence,
            },
        );
    }
    Ok(records)
}
