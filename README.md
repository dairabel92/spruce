# SPrUCE  
**SPrUCE: Sigmoid Pi Requiring Ultraconserved Elements**  


SPrUCE is a tool for estimating nucleotide diversity (π) from Ultraconserved Element (UCE) alignments. It processes MAFFT-aligned UCE loci and models diversity across flanking regions using a Gompertz function, returning per-position π estimates and a final θ (genome-wide π) value.

## Installation

Recommended to create conda environment:

```bash
conda create -n spruce_env python=3.10 -y
conda activate spruce_env
```

Install directly from Github: 

```bash
pip install git+https://github.com/dairabel92/spruce.git
```

After installation, verify:

spruce --help


SPrUCE integrates directly with the standard PHYLUCE pipeline (Faircloth Lab). The recommended input is the output of:

```bash
phyluce_align_seqcap_align \
    --input uce_sequences.fasta \
    --output mafft-edge-trimmed \
    --taxa 10 \
    --aligner mafft \
    --cores 12 \
    --output-format fasta \
    --incomplete-matrix
```

For population genomic analyses, edge trimming is recommended, but internal trimming should be avoided.

## Usage

SPrUCE accepts UCE alignments in: FASTA, FA, FAS, NEXUS, PHY, PHYLIP, CLUSTAL, EMBOSS, STOCKHOLM
Note: For PHYLIP alignments containing long taxon names, the relaxed PHYLIP format is recommended.

Required arguments: 

| Flag | Description |
|------|-------------|
| `--alignments` | Directory containing UCE alignment files |
| `--output-file` | Output file name |
| `--method` | `stack` or `concat` |

---

SPrUCE will take UCE loci  alignment as input and has two method options:

### **STACK**
Stacks (aggregates) per-position nucleotide diversity (π) across all UCEs at each distance from the UCE center. A Gompertz decay curve is then fit to these position-specific means.

### **CONCAT**
Retains locus-specific nucleotide diversity (π) values at each position rather than averaging across UCEs. A single global Gompertz curve is fit by minimizing residuals across all locus–position observations simultaneously.

SPrUCE produces:

- **Final θ estimate** (genome-wide π) printed to screen
- **Per-position π estimates**  written to file

---

##  Optional Arguments

Users may manually control key parameters, including flank size (distance from the UCE center), weighting during model fitting, and minimum coverage thresholds. By default, SPrUCE selects an appropriate flank internally and applies coverage-based weighting, but these settings can be overridden if desired.

| Option | Description |
|--------|-------------|
| `--flank` | Flank size on each side of the UCE core (e.g., 400 or 750) |
| `--use-weights` | TRUE/FALSE (default TRUE). Applies coverage-based weighting during Gompertz fitting |
| `--min-bases` | Minimum number of bases required for a position to be included |

---

## Test Example

Example command using the included test dataset:

```bash
spruce \
  --alignments ./tests/data/plectrophenax_hyperboreus_uce \
  --output-file spruce_test_output.csv \
  --method stack
```
 
## Contact
Developed by **Daira Melendez Belardi, Ali Osman Berk Sapci**  
UC San Diego — Mirarab Lab  
GitHub: https://github.com/dairabel92

