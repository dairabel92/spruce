#SPrUCE  
**SPrUCE: Sigmoid Pi Requiring Ultraconserved Elements**  
Tool for estimating nucleotide diversity (π) from UCE alignments.

---

SPrUCE is a lightweight and fast command-line tool for estimating nucleotide diversity (π) from Ultraconserved Element (UCE) alignments. It processes MAFFT-aligned UCE loci and models diversity across flanking regions, returning per-position π estimates and a final θ value.

SPrUCE integrates directly with the standard PHYLUCE pipeline (Faircloth Lab). The recommended input is the output of:

phyluce_align_seqcap_align \
    --input uce_sequences.fasta \
    --output mafft-edge-trimmed \
    --taxa 10 \
    --aligner mafft \
    --cores 12 \
    --output-format fasta \
    --incomplete-matrix

For population genomic analyses, edge trimming is recommended, but internal trimming should be avoided.

SPrUCE accepts UCE alignments in:
- FASTA  
- NEXUS  

Required arguments:

| Flag | Description |
|------|-------------|
| `--alignments` | Directory containing UCE alignment files |
| `--output-file` | Output file name |
| `--method` | `stack` or `concat` |

---

##  Options

| Option | Description |
|--------|-------------|
| `--flank` | Flank size on each side of the UCE core (e.g., 400 or 750) |
| `--use-weights` | TRUE/FALSE (default TRUE). Applies weight. |
| `--min-base` | Minimum number of bases required for a position to be included |

---


### **STACK**
Stacks (aggregates) per-position nucleotide diversity (π) across UCEs. Fits the Gompertz decay curve.

### **CONCAT**
Fits a Gompertz model to individual UCE locus, treating loci independently rather than stacking by position.

SPrUCE produces:

- **Final θ estimate**, printed to screen and written to file  
- **Per-position π estimates** across all UCEs  
- Optional logs for alignment quality and skipped positions  

---

## example command

python spruce.py
--alignments ./uce_alignments/
--output-file spruce_output.csv
--method stack


This basic command runs SPrUCE using the STACK method and automatically determines: an internal flank threshold based on UCE structure, the minimum number of bases required (--min-base), and default weighting (--use-weights TRUE)

##contact
Developed by **Daira Melendez Belardi, Ali Osman Berk Sapci**  
UC San Diego — Mirarab Lab  
GitHub: https://github.com/dairabel92

