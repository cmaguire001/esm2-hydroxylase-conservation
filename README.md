# esm2-hydroxylase-conservation
Comparative analysis of aromatic amino acid hydroxylase conservation using Meta's ESM-2 protein language model. Validates 2011 manual alignment findings (T-Coffee, MEME, PHYLIP) against transformer-based attention weights. Python · Jupyter · Biopython · PyTorch · UniProt API · Bioinformatics · Computational Biology · Protein Language Models

# Aromatic Amino Acid Hydroxylase Analysis (2011 Replication)

Recreating a 2011 bioinformatics analysis of conserved residues in aromatic amino acid hydroxylases using modern tools in GitHub Codespaces.

## Background

This analysis examines 4 subfamilies of aromatic amino acid hydroxylases:
- **PAH**: Phenylalanine hydroxylase
- **TH**: Tyrosine hydroxylase  
- **TPH**: Tryptophan hydroxylase
- **YOMO**: Tyrosine monooxygenase

### Ground Truth (from 2011 poster)

**11 Invariant Residues:**
- S-349, P-281, G-289, G-344, H-285, H-290, E-330, P-225, R-270, D-282, E-353

**Functional Annotations:**
- **Iron binding triad**: His-285, His-290, Glu-330
- **Salt bridge**: Arg-270, Asp-282

## Quick Start in GitHub Codespaces

### Step 1: Create Repository

1. Go to GitHub and create a new repository: `aromatic-hydroxylases-analysis`
2. Make it public or private (your choice)
3. Initialize with a README (or push this one)
4. Click the green **Code** button → **Codespaces** tab → **Create codespace on main**

### Step 2: Set Up Directory Structure

Once your Codespace opens, run:

```bash
# Create directory structure
mkdir -p data/{raw,processed,results}
mkdir -p scripts
mkdir -p notebooks
mkdir -p docs

# Move the fetch script to scripts folder (if not already there)
mv 01_fetch_uniprot.py scripts/ 2>/dev/null || true
```

### Step 3: Install Dependencies

```bash
# Update package manager
sudo apt-get update

# Install Clustal Omega (for sequence alignment)
sudo apt-get install -y clustalo

# Install Python dependencies
pip install -r requirements.txt
```

### Step 4: Run Data Retrieval

```bash
# Navigate to scripts directory
cd scripts

# Run the UniProt fetcher
python 01_fetch_uniprot.py
```

**Expected output:**
- `data/raw/PAH_sequences.fasta` - Phenylalanine hydroxylase sequences
- `data/raw/TH_sequences.fasta` - Tyrosine hydroxylase sequences
- `data/raw/TPH_sequences.fasta` - Tryptophan hydroxylase sequences
- `data/raw/YOMO_sequences.fasta` - Tyrosine monooxygenase sequences
- `data/raw/all_hydroxylases.fasta` - Combined dataset
- `data/raw/*_metadata.json` - Sequence metadata

### Step 5: Verify Data

```bash
# Check how many sequences were retrieved
ls -lh ../data/raw/

# Quick look at the sequences
head -n 20 ../data/raw/all_hydroxylases.fasta

# Count sequences per subfamily
grep -c "^>" ../data/raw/PAH_sequences.fasta
grep -c "^>" ../data/raw/TH_sequences.fasta
grep -c "^>" ../data/raw/TPH_sequences.fasta
grep -c "^>" ../data/raw/YOMO_sequences.fasta
```

## Project Structure

```
aromatic-hydroxylases-analysis/
├── data/
│   ├── raw/                          # Raw UniProt sequences
│   │   ├── PAH_sequences.fasta
│   │   ├── TH_sequences.fasta
│   │   ├── TPH_sequences.fasta
│   │   ├── YOMO_sequences.fasta
│   │   ├── all_hydroxylases.fasta   # Combined dataset
│   │   └── *_metadata.json          # Sequence metadata
│   ├── processed/                    # Aligned sequences
│   └── results/                      # Analysis outputs
├── scripts/
│   ├── 01_fetch_uniprot.py          # ✓ Phase 1: Data retrieval
│   ├── 02_align_sequences.py        # Phase 2: Multiple sequence alignment
│   └── 03_analyze_conservation.py   # Phase 3: Conservation analysis
├── notebooks/                        # Jupyter notebooks for exploration
├── docs/                             # Documentation
├── requirements.txt                  # Python dependencies
└── README.md                         # This file
```

## Analysis Pipeline

### Phase 1: Data Retrieval (Day 1) ✓
**Status:** Complete  
**Output:** ~20-30 representative sequences from UniProt

### Phase 2: Sequence Alignment (Day 2)
**Goal:** Align sequences using Clustal Omega  
**Output:** Multiple sequence alignment (MSA) in FASTA/Stockholm format

### Phase 3: Conservation Analysis (Day 3)
**Goal:** Identify conserved residues and validate against ground truth  
**Output:** Conservation scores, residue annotations, validation report

## Ground Truth Reference

### Expected Invariant Residues (100% conservation)
Position numbering based on original 2011 analysis:

| Residue | Position | Function |
|---------|----------|----------|
| His-285 | 285 | Iron binding |
| His-290 | 290 | Iron binding |
| Glu-330 | 330 | Iron binding |
| Arg-270 | 270 | Salt bridge |
| Asp-282 | 282 | Salt bridge |
| Ser-349 | 349 | Invariant |
| Pro-281 | 281 | Invariant |
| Gly-289 | 289 | Invariant |
| Gly-344 | 344 | Invariant |
| Pro-225 | 225 | Invariant |
| Glu-353 | 353 | Invariant |

**Note:** Position numbers are from the original 2011 analysis. After alignment, you'll need to map these to alignment columns.

## Next Steps

After completing Phase 1:

1. **Inspect the sequences** - Open FASTA files to verify quality
2. **Check metadata** - Review organism diversity and sequence lengths
3. **Prepare for alignment** - Ensure all sequences are present
4. **Move to Phase 2** - Run multiple sequence alignment

## Troubleshooting

### UniProt API not responding
- Check internet connection
- Verify UniProt REST API is accessible: https://rest.uniprot.org/
- Try reducing `max_seqs` in the script

### No sequences found for YOMO
- This subfamily might be sparse in Swiss-Prot
- Try adjusting the query in `SUBFAMILIES['YOMO']['query']`
- Consider using TrEMBL (remove `reviewed:true`) if needed

### Clustal Omega not found
```bash
sudo apt-get update
sudo apt-get install -y clustalo
clustalo --version  # Should show version info
```

## References

- UniProt REST API: https://www.uniprot.org/help/api
- Clustal Omega: http://www.clustal.org/omega/
- Biopython AlignIO: https://biopython.org/wiki/AlignIO

## Author Notes

This replicates the 2011 poster analysis using a representative subset (~20-30 sequences) instead of all 108 original sequences. The reduced dataset keeps Codespaces runtime manageable while maintaining statistical validity for conservation analysis.

---

**Last Updated:** March 2026  
**Original Analysis:** 2011
