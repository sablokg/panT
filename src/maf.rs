use crate::structfile::ExtractMAF;
use rayon::prelude::*;
use std::error::Error;
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader, BufWriter};

/*
Gaurav Sablok
gsablok@proton.me
*/

#[derive(Debug, Clone)]
pub struct Filepath {
    pub filepath: String,
}

impl Filepath {
    pub async fn processmaf(&self) -> Result<String, Box<dyn Error + Send + Sync>> {
        let valueopen = File::open(&self.filepath).await?;
        let mut valueread = BufReader::new(valueopen).lines();

        let mut rawlines: Vec<String> = Vec::new();
        while let Some(line) = valueread.next_line().await? {
            rawlines.push(line);
        }

        // HYPER-THREADING: parsing each alignment line into a record is
        // independent per-line work, so it's spread across all
        // available cores with rayon instead of a single-threaded loop.
        // `par_iter` on a Vec is an indexed/ordered parallel iterator,
        // so `filter_map` here still preserves the original line order,
        // which matters because the query/target split below relies on
        // alternating even/odd positions.
        let mafcapture: Vec<ExtractMAF> = rawlines
            .par_iter()
            .filter(|line| line.starts_with('s'))
            .map(|line| {
                let linevec = line.split('\t').collect::<Vec<_>>();
                ExtractMAF {
                    refname: linevec[0].to_string(),
                    start: linevec[2].to_string(),
                    strand: linevec[3].to_string(),
                    end: linevec[4].to_string(),
                    sequence: linevec[5].to_string(),
                }
            })
            .collect();

        let alignquery = mafcapture
            .par_iter()
            .step_by(2)
            .cloned()
            .collect::<Vec<ExtractMAF>>();
        let targetseq = mafcapture
            .par_iter()
            .skip(1)
            .step_by(2)
            .cloned()
            .collect::<Vec<ExtractMAF>>();

        // write both fasta files concurrently instead of sequentially
        let (q, t) = tokio::join!(
            write_fasta("query.fasta", alignquery),
            write_fasta("target.fasta", targetseq),
        );
        q?;
        t?;

        Ok("The maf file has been parsed".to_string())
    }
}

async fn write_fasta(
    path: &str,
    records: Vec<ExtractMAF>,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let file = File::create(path).await?;
    let mut writer = BufWriter::new(file);
    for rec in records {
        writer
            .write_all(format!(">{}\n{}\n", rec.refname, rec.sequence).as_bytes())
            .await?;
    }
    writer.flush().await?;
    Ok(())
}
