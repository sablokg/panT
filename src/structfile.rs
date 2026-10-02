#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct PAFAlign {
    pub queryname: String,
    pub querystart: String,
    pub queryend: String,
    pub refname: String,
    pub refstart: String,
    pub refend: String,
}

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct SamAlign {
    pub refname: String,
    pub queryname: String,
    pub querystart: usize,
    pub queryend: usize,
    pub refstart: usize,
    pub refend: usize,
    pub queryseq: String,
}

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct ExtractMAF {
    pub refname: String,
    pub start: String,
    pub strand: String,
    pub end: String,
    pub sequence: String,
}

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct PAFStrand {
    pub queryname: String,
    pub querystart: String,
    pub queryend: String,
    pub querystrand: String,
    pub refname: String,
    pub refstart: String,
    pub refend: String,
    pub refstrand: String,
}
