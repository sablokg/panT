use clap::{Parser, Subcommand};
#[derive(Debug, Parser)]
#[command(
    name = "pantranscript",
    version = "1.0",
    about = "alignment extraction for pantranscriptomics
       ************************************************
       Author Gaurav Sablok
       Email gsablok@proton.me
      ************************************************"
)]
pub struct CommandParse {
    /// subcommands for the specific actions
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// SAM alignment
    SAMAlignment {
        /// path to the aligned sam
        pathname: String,
        /// path to the referencefasta
        fasta: String,
    },
    /// MAF aligned file
    MAFAlignment {
        /// path to the alignment
        alignment: String,
    },
    /// PAF alignment filter
    Pafalignment {
        /// paf alignment
        pafalignment: String,
        /// filter residue number
        filterresiduematch: String,
    },
    /// PAF alignment slicer
    Pafalign {
        /// paf aligned file
        pafaligned: String,
        /// reference fasta file
        reference: String,
        /// query fasta file
        query: String,
    },
    /// PAF filter strand
    PafStrand {
        /// path to the aligned file
        alignedfile: String,
        /// strand type
        strand: Option<String>,
    },
    /// PAF filter strand and aligned bases
    PafStrandbases {
        /// path to the aligned file
        alignedfile: String,
        /// strand type
        strand: String,
        /// number of bases
        bases: String,
    },
}
