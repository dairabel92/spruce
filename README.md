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

```bash
spruce --help
```

SPrUCE can integrate directly with the standard PHYLUCE pipeline (Faircloth Lab). The recommended input is the output of this step:

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

UCE alignments can be generated from: 
- Targeted UCE sequence capture data (e.g., biological samples enriched with UCE probes);
- Whole-genome assemblies (e.g., WGS data where UCE loci are extracted *in silico*).

For details on how to extract and align UCE loci from raw reads or genome assemblies, see the [PHYLUCE tutorials](https://phyluce.readthedocs.io/en/latest/tutorials/index.html).

## Usage

SPrUCE accepts UCE alignments in: FASTA, FA, FAS, NEXUS, PHY, PHYLIP, CLUSTAL, EMBOSS, STOCKHOLM

Note: For PHYLIP alignments containing long taxon names, the relaxed PHYLIP format (`phylip-relaxed`) is recommended to ensure compatibility with BioPython/AlignIO parsing.

### Required arguments: 

| Flag | Description |
|------|-------------|
| `--alignments` | Directory containing UCE alignment files |
| `--output-file` | Output file name for csv file of all sites|
| `--output-estimate` | Output file name for final estimate of theta and other gompertz parameters |
| `--method` | `stack` (default) or `concat` |

---

SPrUCE will take UCE loci alignment as input and has two method options:

### **STACK** (Default)
Stacks (aggregates) per-position nucleotide diversity (π) across all UCEs at each distance from the UCE center. A Gompertz decay curve is then fit to these position-specific means.

### **CONCAT**
Retains locus-specific nucleotide diversity (π) values at each position rather than averaging across UCEs. A Gompertz curve is fit by minimizing residuals across all locus–position observations simultaneously.

> **Please note:** In the current implementation, increasing the number of cores does **not** speed up `CONCAT` mode.  
>
> We strongly recommend submitting `CONCAT` runs as a cluster job (e.g., via SLURM) rather than running interactively, as this mode can take substantially longer to complete depending on dataset size.


### **SPrUCE PRODUCES:**

- **Final θ estimate** (genome-wide π) printed to screen and saved to `--output-estimate` file
- **Per-position π estimates**  written to `--output-file` file 

---

##  Optional Arguments

Users may manually control key parameters, including flank size (distance from the UCE center), weighting during model fitting, and minimum coverage thresholds. By default, SPrUCE selects an appropriate flank internally and applies coverage-based weighting, but these settings can be overridden if desired.

| Option | Description                                                                                                                                                                                                                                    |
|--------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--flank` | Flank size on each side of the UCE core (e.g., 400 or 750). Positions beyond this distance from the center are ignored. Default: automatically determined from the data.                                  |
| `--use-weights` | TRUE/FALSE. Default: TRUE. Applies sample-size–based weighting during Gompertz fitting to down-weight positions with low sample size (e.g., gappy or missing sites). |                         
| `--min-bases` | Minimum number of bases required for a position to be included; if a site is gappy has fewer than this provided number of non-gap sequences, we ignore it. Default: automatically determined based on coverage. |

---

## Test Example

The test dataset consists of UCE alignments from McKay's Bunting (*Plectrophenax hyperboreus*) published in Winker et al. (2018) (see citation below).

Clone the GitHub repository to access test dataset:

```bash
git clone https://github.com/dairabel92/spruce.git
```

Example command using the included test dataset:

```bash
spruce \
  --alignments ./tests/data/plectrophenax_hyperboreus_uce \
  --output-estimate estimates.csv \ 
  --output-file spruce_test_output.csv \
  --method stack
```
 
## Contact
Developed by **Daira Melendez Belardi, Ali Osman Berk Sapci**  
UC San Diego — Mirarab Lab  
GitHub: https://github.com/dairabel92

---

## Data Source

The included test dataset is derived from the UCE data published in:

Winker, K., Glenn, T. C., Faircloth, B. C., *et al.* (2018). *Ultraconserved elements (UCEs) illuminate the population genomics of a recent, high-latitude avian speciation event*. *PeerJ* **6**, e5735. doi:10.7717/peerj.5735

Please cite the original publication if using these data.

