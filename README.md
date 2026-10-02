# pantranscript

- a uniform software for all the pantranscriptomics
- gives you the aligned region from several alignment format such as SAM, PAF, MAF for aligned transcriptome to the reference transcriptome
- build for my entire pipeline for the pangenomes https://codeberg.org/gsablok/graphpan and https://codeberg.org/gsablok/panseq

```
cargo build
```

```
                                                                    
                                                        _|_|_|_|_|  
                       _|_|_|       _|_|_|   _|_|_|         _|      
                       _|    _|   _|    _|   _|    _|       _|      
                       _|    _|   _|    _|   _|    _|       _|      
                       _|_|_|       _|_|_|   _|    _|       _|      
                       _|                                           
                       _|                                           
alignment extraction for pantranscriptomics
       ************************************************
       Author Gaurav Sablok
       Email gsablok@proton.me
      ************************************************

Usage: pant <COMMAND>

Commands:
  sam-alignment    SAM alignment
  maf-alignment    MAF aligned file
  pafalignment     PAF alignment filter
  pafalign         PAF alignment slicer
  paf-strand       PAF filter strand
  paf-strandbases  PAF filter strand and aligned bases
  help             Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

All of the parsing/slicing/filtering work in every command below is spread
across all available CPU cores with [rayon](https://docs.rs/rayon), so
larger alignment files scale with the number of cores on your machine.

---

## `sam-alignment`

Extracts, from a SAM alignment, the aligned reference region and the
aligned query sequence for every record, and writes them out as two
matching fasta files.

```
pant sam-alignment <PATHNAME> <FASTA>
```

| Argument | Meaning |
|---|---|
| `PATHNAME` | path to the aligned SAM file |
| `FASTA` | path to the reference fasta the SAM file was aligned against |

**Produces:** `reference.fasta`, `query.fasta`

```bash
pant sam-alignment test/aln.sam test/sample.fa
```

For each non-header SAM line, `pantranscript` computes the reference
window (`POS - 1 .. POS - 1 + len(SEQ)`) from the reference fasta and
pairs it with the query sequence from column 10 of the SAM record.

---

## `maf-alignment`

Splits a MAF (Multiple Alignment Format) file into its two alternating
alignment blocks — query and target — and writes each out as a fasta
file.

```
pant maf-alignment <ALIGNMENT>
```

| Argument | Meaning |
|---|---|
| `ALIGNMENT` | path to the `.maf` alignment file |

**Produces:** `query.fasta`, `target.fasta`

```bash
pant maf-alignment test/aln.maf
```

Every line beginning with `s` in the MAF file is one aligned sequence
block; `pantranscript` takes every other block as the query and the
remaining ones as the target (the standard two-line-per-block MAF
layout).

---

## `pafalignment`

Filters a PAF alignment file down to records whose residue-match count
(PAF column 10) is less than or equal to a given threshold, and writes
the matching records back out in PAF format.

```
pant pafalignment <PAFALIGNMENT> <FILTERRESIDUEMATCH>
```

| Argument | Meaning |
|---|---|
| `PAFALIGNMENT` | path to the `.paf` alignment file |
| `FILTERRESIDUEMATCH` | keep records with `residues <= this number` |

**Produces:** `filtered-paf.txt`

```bash
pant pafalignment test/K-3138.paf 3000
```

---

## `pafalign`

Slices out the aligned region of both the reference and query fasta
files for every record in a PAF file, using the PAF start/end
coordinates (columns 3/4 for the query, columns 8/9 for the reference).

```
pant pafalign <PAFALIGNED> <REFERENCE> <QUERY>
```

| Argument | Meaning |
|---|---|
| `PAFALIGNED` | path to the `.paf` alignment file |
| `REFERENCE` | reference fasta the alignment was made against |
| `QUERY` | query fasta the alignment was made from |

**Produces:** `reference.fasta`, `query.fasta`

```bash
pant pafalign test/aln.paf test/shorties.fa test/sample.fa
```

The reference and query fasta files are read concurrently, and every
record's aligned slice is looked up and cut out in parallel.

---

## `paf-strand`

Splits a PAF file into plus-strand and minus-strand alignments (PAF
column 5, the single strand field for the whole alignment).

```
pant paf-strand <ALIGNEDFILE> [STRAND]
```

| Argument | Meaning |
|---|---|
| `ALIGNEDFILE` | path to the `.paf` alignment file |
| `STRAND` | optional: `plus` to only compute the plus-strand file, `neg` to only compute the minus-strand file. Omit to skip both. |

**Produces:** `plustrand.txt` (when `STRAND=plus`), `negativestrand.txt`
(when `STRAND=neg`)

```bash
pant paf-strand test/K-3138.paf plus
pant paf-strand test/K-3138.paf neg
```

---

## `paf-strandbases`

Filters a PAF file by strand *and* a residue-match threshold at the
same time, writing the matching records back out in PAF format.

```
pant paf-strandbases <ALIGNEDFILE> <STRAND> <BASES>
```

| Argument | Meaning |
|---|---|
| `ALIGNEDFILE` | path to the `.paf` alignment file |
| `STRAND` | strand to keep, e.g. `+` or `-` |
| `BASES` | keep records with `residues <= this number` |

**Produces:** `filtered-paf.txt`

```bash
pant paf-strandbases test/K-3138.paf + 3000
```

---

Gaurav Sablok \
codeprog@icloud.com
