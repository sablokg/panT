use crate::structfile::PAFAlign;
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::io::Write;
use std::io::{self, BufRead};
use std::path::Path;

/*
Gaurav Sablok
gsablok@proton.me
*/

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Filename {
    pub aligned: String,
    pub reference: String,
    pub queryname: String,
}

impl Filename {
    pub fn parsepaf(&self) -> Result<String, Box<dyn Error>> {
        let alignedfile = self.aligned.clone();

        // BUG FIX: `read_fasta` used to be an `async fn` wrapped in
        // `#[tokio::main]`. Calling it here (itself invoked from inside
        // main's own `#[tokio::main]` runtime) panicked with "Cannot
        // start a runtime from within a runtime" every single time this
        // subcommand ran. `read_fasta` never actually awaited anything,
        // so it is now a plain synchronous function.
        //
        // HYPER-THREADING: the reference and query fasta files are
        // completely independent, so read them on two threads at once
        // with rayon::join instead of one after another.
        let (referencefasta, queryfasta) = rayon::join(
            || read_fasta(self.reference.clone()),
            || read_fasta(self.queryname.clone()),
        );
        let referencefasta = referencefasta?;
        let queryfasta = queryfasta?;

        let mut pafalign: Vec<PAFAlign> = Vec::new();
        let alignopen = File::open(alignedfile).expect("file not present");
        let alignedread = BufReader::new(alignopen);
        for i in alignedread.lines() {
            let line = i.expect("line not present");
            let linevec = line.split("\t").collect::<Vec<_>>();
            pafalign.push(PAFAlign {
                queryname: linevec[0].to_string(),
                querystart: linevec[2].to_string(),
                queryend: linevec[3].to_string(),
                refname: linevec[5].to_string(),
                refstart: linevec[7].to_string(),
                refend: linevec[8].to_string(),
            })
        }

        // BUG FIX: the original code re-scanned the *entire* fasta
        // HashMap with a nested loop for every single PAF record
        // (an O(n * m) linear search) instead of doing a direct
        // key lookup, since `queryfasta`/`referencefasta` are already
        // keyed HashMaps.
        //
        // HYPER-THREADING: slicing out each aligned region is
        // independent per-record work, so it is done with rayon's
        // parallel iterator across all available cores. `par_iter`
        // on a Vec preserves order, so the output files come out in
        // the same order as before.
        let query: Vec<(String, String)> = pafalign
            .par_iter()
            .filter_map(|i| {
                queryfasta.get(&i.queryname).and_then(|value| {
                    let start = i.querystart.parse::<usize>().ok()?;
                    let end = i.queryend.parse::<usize>().ok()?;
                    value
                        .sequence
                        .get(start..end)
                        .map(|s| (i.queryname.clone(), s.to_string()))
                })
            })
            .collect();

        let reference: Vec<(String, String)> = pafalign
            .par_iter()
            .filter_map(|i| {
                referencefasta.get(&i.refname).and_then(|value| {
                    let start = i.refstart.parse::<usize>().ok()?;
                    let end = i.refend.parse::<usize>().ok()?;
                    value
                        .sequence
                        .get(start..end)
                        .map(|s| (i.refname.clone(), s.to_string()))
                })
            })
            .collect();

        // Write both output files concurrently instead of sequentially.
        rayon::join(
            || {
                let mut filequery = File::create("query.fasta").expect("file not present");
                for i in query.iter() {
                    writeln!(filequery, ">{}\n{}", i.0, i.1).expect("file not present");
                }
            },
            || {
                let mut fileref = File::create("reference.fasta").expect("file not found");
                for i in reference.iter() {
                    writeln!(fileref, ">{}\n{}", i.0, i.1).expect("file not present");
                }
            },
        );

        Ok("PAF has been parsed".to_string())
    }
}

#[derive(Debug)]
pub struct FastaRecord {
    pub id: String,
    pub sequence: String,
}

/// Plain synchronous fasta reader. Previously this was declared as
/// `async fn` and annotated with `#[tokio::main]`, which made it spin up
/// its own Tokio runtime every time it was called. Since it was always
/// called from *inside* `main`'s own `#[tokio::main]` runtime, that
/// caused a guaranteed panic ("Cannot start a runtime from within a
/// runtime") the moment `pafalign` was run. The body never awaited
/// anything anyway, so it is simply a regular function now.
pub fn read_fasta<P: AsRef<Path>>(path: P) -> io::Result<HashMap<String, FastaRecord>> {
    let file = File::open(path)?;
    let reader = io::BufReader::new(file);
    let mut records = HashMap::new();
    let mut current_id = String::new();
    let mut current_sequence = String::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_id.is_empty() {
                records.insert(
                    current_id.clone(),
                    FastaRecord {
                        id: current_id.clone(),
                        sequence: current_sequence.clone(),
                    },
                );
                current_sequence.clear();
            }
            current_id = line[1..].to_string();
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
