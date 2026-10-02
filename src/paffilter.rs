use rayon::prelude::*;
use std::error::Error;
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader};

/*
Gaurav Sablok
gsablok@proton.me
 */

#[derive(Debug, Clone)]
pub struct PAFfile {
    pub filename: String,
}

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct PAFFilter {
    qname: String,
    qseqlength: String,
    qstart: String,
    qend: String,
    qstrand: String,
    target: String,
    targetseqlength: String,
    targetstart: String,
    targetend: String,
    residues: String,
}

impl PAFfile {
    pub async fn filterpaf(&self, residuecount: &str) -> Result<String, Box<dyn Error>> {
        let fileopen = File::open(self.filename.clone()).await?;
        let fileread = BufReader::new(fileopen);
        let mut lines = fileread.lines();

        let mut rawlines: Vec<String> = Vec::new();
        while let Some(line) = lines.next_line().await? {
            rawlines.push(line);
        }

        // HYPER-THREADING: parsing each line into a record is
        // independent per-line work, so it is spread across all
        // available cores with rayon instead of a single-threaded loop.
        let extractedpaf: Vec<PAFFilter> = rawlines
            .par_iter()
            .map(|line| {
                let linevec = line.split('\t').collect::<Vec<_>>();
                PAFFilter {
                    qname: linevec[0].to_string(),
                    qseqlength: linevec[1].to_string(),
                    qstart: linevec[2].to_string(),
                    qend: linevec[3].to_string(),
                    qstrand: linevec[4].to_string(),
                    target: linevec[5].to_string(),
                    targetseqlength: linevec[6].to_string(),
                    targetstart: linevec[7].to_string(),
                    targetend: linevec[8].to_string(),
                    residues: linevec[9].to_string(),
                }
            })
            .collect();

        let residuecount: usize = residuecount.parse()?;

        // HYPER-THREADING: the filter predicate parses a number on
        // every record, so run it in parallel too.
        let filteredpaf: Vec<&PAFFilter> = extractedpaf
            .par_iter()
            .filter(|i| {
                i.residues
                    .parse::<usize>()
                    .map(|r| r <= residuecount)
                    .unwrap_or(false)
            })
            .collect();

        let mut filewrite = File::create("filtered-paf.txt").await?;
        for i in filteredpaf.iter() {
            let line = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                i.qname,
                i.qseqlength,
                i.qstart,
                i.qend,
                i.qstrand,
                i.target,
                i.targetseqlength,
                i.targetstart,
                i.targetend,
                i.residues
            );
            filewrite.write_all(line.as_bytes()).await?;
        }

        Ok("The paf has been filtered".to_string())
    }
}
