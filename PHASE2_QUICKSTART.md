# Phase 2: Multiple Sequence Alignment - Quick Start

## What Phase 2 Does

Aligns your 30 sequences using **Clustal Omega** to:
- Identify conserved positions across all subfamilies
- Create a multiple sequence alignment (MSA)
- Calculate conservation scores for each column
- Prepare data for Phase 3 validation

## Prerequisites

Make sure Clustal Omega is installed:

```bash
clustalo --version
```

If not installed:
```bash
sudo apt-get update
sudo apt-get install -y clustalo
```

## Run Phase 2

From your repository root:

```bash
cd /workspaces/esm2-hydroxylase-conservation

# Copy the script to scripts directory
cp 02_align_sequences.py scripts/

# Run alignment
cd scripts
python 02_align_sequences.py
```

## Expected Output

The script will:
1. Load 30 sequences from `data/raw/all_hydroxylases.fasta`
2. Run Clustal Omega (~30-60 seconds)
3. Create alignment files in `data/processed/`
4. Calculate conservation statistics
5. Identify fully conserved positions

**Output files:**
- `data/processed/alignment.fasta` - Aligned sequences (FASTA format)
- `data/processed/alignment.clustal` - Human-readable alignment
- `data/processed/alignment_stats.json` - Statistics and conservation scores

## What to Expect

**Alignment length:** ~520-550 columns (longer than average sequence due to gaps)

**Conservation:** Should find several fully conserved positions including:
- Iron binding residues (His, His, Glu)
- Salt bridge pair (Arg, Asp)
- Structural residues (Pro, Gly, Ser, etc.)

**Runtime:** 30-60 seconds for 30 sequences

## Verification

After running:

```bash
# Check outputs were created
ls -lh ../data/processed/

# Quick peek at alignment
head -n 20 ../data/processed/alignment.fasta

# Check statistics
cat ../data/processed/alignment_stats.json | grep -A 5 "fully_conserved"
```

## Troubleshooting

**"clustalo: command not found"**
```bash
sudo apt-get install clustalo
```

**"Input file not found"**
```bash
# Verify Phase 1 data exists
ls -lh data/raw/all_hydroxylases.fasta
```

**Alignment takes too long (>2 minutes)**
- This is unusual for 30 sequences
- Check if Codespace has sufficient resources
- Try reducing thread count in script (change `--threads 2` to `--threads 1`)

## Next Step

Once Phase 2 completes successfully:
- **Phase 3**: Conservation analysis & validation
- Map conserved positions to your ground truth (11 invariant residues)
- Generate validation report

---

**Ready to run?** Execute the commands above and Phase 2 should complete in ~1 minute!
