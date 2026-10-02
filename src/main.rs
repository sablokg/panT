mod args;
use crate::args::CommandParse;
use crate::args::Commands;
use clap::Parser;
mod structfile;
use crate::sam::SAMFile;
mod maf;
mod sam;
use self::paffilter::PAFfile;
use crate::maf::Filepath;
use vornix_banner::{Banner, BuiltinFont, Style, rgb};
mod paf;
mod paffilter;
use crate::paf::Filename;
mod strand;
use crate::strand::PAFAlignedFile;
mod pafstrandbases;
use crate::pafstrandbases::PAFfileStrand;

/*
Gaurav Sablok
gsablok@proton.me
*/

#[tokio::main]
async fn main() {
    let style = Style::new().fg(rgb(255, 100, 20)).bold();

    let mut banner = Banner::new("panT")
        .with_builtin_font(BuiltinFont::Block)
        .with_style(style)
        .centered(true);

    banner.display().unwrap();
    let argparse = CommandParse::parse();
    match &argparse.command {
        Commands::SAMAlignment { pathname, fasta } => {
            let command = SAMFile {
                alignedsame: pathname.clone(),
                referencefasta: fasta.clone(),
            };
            let commandnew = command.parsesam().await.unwrap();
            println!("The sam file has been parsed:{}", commandnew);
        }
        Commands::MAFAlignment { alignment } => {
            let value = Filepath {
                filepath: alignment.clone(),
            };
            let command = value.processmaf().await.unwrap();
            println!("The command has finished: {}", command);
        }
        Commands::Pafalignment {
            pafalignment,
            filterresiduematch,
        } => {
            let filename = PAFfile {
                filename: pafalignment.clone(),
            };

            let command = filename.filterpaf(filterresiduematch).await.unwrap();
            println!("The command has finished:{}", command);
        }
        Commands::Pafalign {
            pafaligned,
            reference,
            query,
        } => {
            let command = Filename {
                aligned: pafaligned.clone(),
                reference: reference.clone(),
                queryname: query.clone(),
            };
            let commandnew = command.parsepaf().unwrap();
            println!("The command has been completed: {}", commandnew);
        }
        Commands::PafStrand {
            alignedfile,
            strand,
        } => {
            let command = PAFAlignedFile {
                filename: alignedfile.clone(),
            };
            let runcommand = command.strandaware(strand.clone()).unwrap();
            println!("The command has finished:{}", runcommand);
        }
        Commands::PafStrandbases {
            alignedfile,
            strand,
            bases,
        } => {
            let file = PAFfileStrand {
                filename: alignedfile.clone(),
            };
            let command = file.filterpafstrandbases(bases, strand).await.unwrap();
            println!("The command has finished:{}", command);
        }
    }
}
