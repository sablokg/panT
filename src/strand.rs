use crate::structfile::PAFStrand;
use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};

/*
Gaurav Sablok
gsablok@proton.me
 */

#[derive(Debug, Clone)]
pub struct PAFAlignedFile {
    pub filename: String,
}

impl PAFAlignedFile {
    pub fn strandaware(&self, strand: Option<String>) -> Result<String, Box<dyn Error>> {
        let fileopen = File::open(self.filename.clone()).expect("file not present");
        let fileread = BufReader::new(fileopen);

        let lines: Vec<String> = fileread
            .lines()
            .collect::<Result<Vec<String>, _>>()
            .expect("line not found");

        // BUG FIX: the column indices used here did not match the PAF
        // format at all. Standard PAF columns are:
        //   0 qname, 1 qlen, 2 qstart, 3 qend, 4 strand,
        //   5 tname, 6 tlen, 7 tstart, 8 tend, 9 matches, ...
        // The old code read querystart from column 1 (actually qlen),
        // queryend from column 2 (actually qstart), a "querystrand"
        // from column 3 (actually qend), refname from column 4
        // (actually the strand column), refstart/refend from columns
        // 5/6 (actually tname/tlen), and skipped column 7 entirely to
        // read a "refstrand" from column 8 (actually tend, a number,
        // never "+"/"-"). Because of that last one, `i.refstrand == "+"`
        // could never be true, so every "plus"/"neg" filter silently
        // produced empty output. PAF also only carries a single strand
        // field for the whole alignment (column 4), not separate query
        // and reference strands, so both fields below are filled from
        // that one column.
        //
        // HYPER-THREADING: parsing each line into a record is
        // independent work, so it's done in parallel with rayon.
        let resultvec: Vec<PAFStrand> = lines
            .par_iter()
            .map(|line| {
                let linevec = line.split('\t').collect::<Vec<_>>();
                PAFStrand {
                    queryname: linevec[0].to_string(),
                    querystart: linevec[2].to_string(),
                    queryend: linevec[3].to_string(),
                    querystrand: linevec[4].to_string(),
                    refname: linevec[5].to_string(),
                    refstart: linevec[7].to_string(),
                    refend: linevec[8].to_string(),
                    refstrand: linevec[4].to_string(),
                }
            })
            .collect();

        let want_plus = strand.as_deref() == Some("plus");
        let want_neg = strand.as_deref() == Some("neg");

        // HYPER-THREADING: filtering and writing out the plus-strand and
        // negative-strand results are independent of each other, so run
        // them concurrently, and filter each set in parallel with rayon.
        rayon::join(
            || {
                let plusstrand: Vec<&PAFStrand> = if want_plus {
                    resultvec.par_iter().filter(|i| i.refstrand == "+").collect()
                } else {
                    Vec::new()
                };
                let mut filewrite = File::create("plustrand.txt").expect("file not present");
                for i in plusstrand.iter() {
                    writeln!(
                        filewrite,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        i.queryname,
                        i.querystart,
                        i.queryend,
                        i.querystrand,
                        i.refname,
                        i.refstart,
                        i.refend,
                        i.refstrand
                    )
                    .expect("file not found");
                }
            },
            || {
                let negativestrand: Vec<&PAFStrand> = if want_neg {
                    resultvec.par_iter().filter(|i| i.refstrand == "-").collect()
                } else {
                    Vec::new()
                };
                // BUG FIX: this output filename was misspelled
                // "negaitvestrand.txt", which made it easy to miss the
                // file when scripting against this tool's output.
                let mut negativewrite =
                    File::create("negativestrand.txt").expect("file not present");
                for i in negativestrand.iter() {
                    writeln!(
                        negativewrite,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        i.queryname,
                        i.querystart,
                        i.queryend,
                        i.querystrand,
                        i.refname,
                        i.refstart,
                        i.refend,
                        i.refstrand
                    )
                    .expect("file not found");
                }
            },
        );

        Ok("Result have been parsed according to the strand".to_string())
    }
}
