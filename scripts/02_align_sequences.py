#!/usr/bin/env python3
"""
Phase 2: Multiple Sequence Alignment
Align aromatic amino acid hydroxylase sequences using Clustal Omega.

Input: data/raw/all_hydroxylases.fasta (30 sequences)
Output: data/processed/alignment.fasta (aligned sequences)
        data/processed/alignment.clustal (human-readable format)
        data/processed/alignment_stats.json (alignment statistics)
"""

import subprocess
import json
from pathlib import Path
from collections import Counter
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

# Get repository root
SCRIPT_DIR = Path(__file__).parent.resolve()
REPO_ROOT = SCRIPT_DIR.parent if SCRIPT_DIR.name == 'scripts' else SCRIPT_DIR

# Paths
DATA_RAW = REPO_ROOT / "data" / "raw"
DATA_PROCESSED = REPO_ROOT / "data" / "processed"
DATA_PROCESSED.mkdir(parents=True, exist_ok=True)

INPUT_FASTA = DATA_RAW / "all_hydroxylases.fasta"
OUTPUT_FASTA = DATA_PROCESSED / "alignment.fasta"
OUTPUT_CLUSTAL = DATA_PROCESSED / "alignment.clustal"
OUTPUT_STATS = DATA_PROCESSED / "alignment_stats.json"

print("=" * 80)
print("Phase 2: Multiple Sequence Alignment")
print("=" * 80)
print()
print(f"Repository root: {REPO_ROOT}")
print(f"Input: {INPUT_FASTA}")
print(f"Output directory: {DATA_PROCESSED}")
print()

# Verify input exists
if not INPUT_FASTA.exists():
    print(f"❌ ERROR: Input file not found: {INPUT_FASTA}")
    print("Run Phase 1 (01_fetch_uniprot.py) first!")
    exit(1)

# Count input sequences
input_seqs = list(SeqIO.parse(INPUT_FASTA, "fasta"))
print(f"✓ Found {len(input_seqs)} sequences to align")
print()

# Check sequence length distribution before alignment
seq_lengths = [len(seq.seq) for seq in input_seqs]
print("Input sequence statistics:")
print(f"  Shortest: {min(seq_lengths)} aa")
print(f"  Longest: {max(seq_lengths)} aa")
print(f"  Average: {sum(seq_lengths)/len(seq_lengths):.1f} aa")
print()

# Run Clustal Omega
print("=" * 80)
print("Running Clustal Omega...")
print("=" * 80)
print("This may take 30-60 seconds for 30 sequences...")
print()

try:
    # Run clustalo command
    cmd = [
        "clustalo",
        "-i", str(INPUT_FASTA),
        "-o", str(OUTPUT_FASTA),
        "--outfmt", "fasta",
        "--force",  # Overwrite if exists
        "--verbose",  # Show progress
        "--threads", "2"  # Use 2 threads for speed
    ]
    
    print(f"Command: {' '.join(cmd)}")
    print()
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=300  # 5 minute timeout
    )
    
    if result.returncode != 0:
        print(f"❌ Clustal Omega failed!")
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        exit(1)
    
    print("✓ Alignment complete!")
    print()
    
    # Also save in Clustal format for readability
    print("Creating Clustal format output for viewing...")
    cmd_clustal = [
        "clustalo",
        "-i", str(INPUT_FASTA),
        "-o", str(OUTPUT_CLUSTAL),
        "--outfmt", "clustal",
        "--force"
    ]
    subprocess.run(cmd_clustal, check=True, capture_output=True)
    print(f"✓ Saved: {OUTPUT_CLUSTAL}")
    
except subprocess.TimeoutExpired:
    print("❌ Alignment timed out (>5 minutes)")
    exit(1)
except subprocess.CalledProcessError as e:
    print(f"❌ Error running Clustal Omega: {e}")
    exit(1)
except FileNotFoundError:
    print("❌ Clustal Omega not found!")
    print("Install with: sudo apt-get install clustalo")
    exit(1)

# Parse and analyze alignment
print()
print("=" * 80)
print("Analyzing Alignment")
print("=" * 80)
print()

alignment = AlignIO.read(OUTPUT_FASTA, "fasta")
alignment_length = alignment.get_alignment_length()
num_sequences = len(alignment)

print(f"Sequences: {num_sequences}")
print(f"Alignment length: {alignment_length} columns")
print()

# Calculate statistics
total_gaps = sum(str(record.seq).count('-') for record in alignment)
total_positions = num_sequences * alignment_length
gap_percentage = (total_gaps / total_positions) * 100

print(f"Total positions: {total_positions:,}")
print(f"Gap positions: {total_gaps:,}")
print(f"Gap percentage: {gap_percentage:.2f}%")
print()

# Calculate per-column statistics
print("Calculating conservation scores...")

column_stats = []
fully_conserved_positions = []
gap_heavy_columns = []

for i in range(alignment_length):
    column = alignment[:, i]
    
    # Count residues (excluding gaps)
    residues = [res for res in column if res != '-']
    
    if not residues:
        # All gaps
        conservation = 0.0
        consensus = '-'
    else:
        # Most common residue
        counter = Counter(residues)
        consensus, count = counter.most_common(1)[0]
        
        # Conservation score: fraction of non-gap positions with consensus
        conservation = count / len(residues) if residues else 0.0
        
        # Check for full conservation
        if conservation == 1.0 and len(residues) == num_sequences:
            fully_conserved_positions.append({
                'column': i + 1,  # 1-indexed
                'residue': consensus
            })
    
    # Track gap-heavy columns
    gap_count = column.count('-')
    if gap_count / num_sequences > 0.5:
        gap_heavy_columns.append(i + 1)
    
    column_stats.append({
        'column': i + 1,
        'consensus': consensus,
        'conservation': conservation,
        'gaps': gap_count
    })

print(f"✓ Found {len(fully_conserved_positions)} fully conserved positions (100% identity)")
print(f"✓ Found {len(gap_heavy_columns)} gap-heavy columns (>50% gaps)")
print()

# Display first few conserved positions
if fully_conserved_positions:
    print("Preview of fully conserved positions:")
    for pos in fully_conserved_positions[:10]:
        print(f"  Column {pos['column']:4d}: {pos['residue']}")
    if len(fully_conserved_positions) > 10:
        print(f"  ... and {len(fully_conserved_positions) - 10} more")
    print()

# Save statistics
stats = {
    'input_file': str(INPUT_FASTA),
    'output_file': str(OUTPUT_FASTA),
    'num_sequences': num_sequences,
    'alignment_length': alignment_length,
    'total_positions': total_positions,
    'gap_positions': total_gaps,
    'gap_percentage': gap_percentage,
    'fully_conserved_count': len(fully_conserved_positions),
    'fully_conserved_positions': fully_conserved_positions,
    'gap_heavy_columns': gap_heavy_columns,
    'column_statistics': column_stats  # All columns with conservation scores
}

with open(OUTPUT_STATS, 'w') as f:
    json.dump(stats, f, indent=2)

print(f"✓ Saved statistics: {OUTPUT_STATS}")
print()

# Preview alignment
print("=" * 80)
print("Alignment Preview (first 3 sequences, first 80 columns)")
print("=" * 80)
print()

for i, record in enumerate(alignment[:3]):
    # Get subfamily from header
    header_parts = record.id.split('|')
    accession = header_parts[0] if len(header_parts) > 0 else "???"
    subfamily = header_parts[2] if len(header_parts) > 2 else "???"
    organism = header_parts[3] if len(header_parts) > 3 else "???"
    
    # Show first 80 positions
    seq_preview = str(record.seq)[:80]
    
    print(f"{accession:12s} ({subfamily:4s}) {organism[:20]:20s}")
    print(f"{'':45s}{seq_preview}")
    print()

print(f"(Full alignment saved to {OUTPUT_FASTA})")
print()

# Summary
print("=" * 80)
print("Phase 2 Complete!")
print("=" * 80)
print()
print("✓ Multiple sequence alignment created")
print(f"✓ Alignment length: {alignment_length} columns")
print(f"✓ Fully conserved positions: {len(fully_conserved_positions)}")
print()
print("Output files:")
print(f"  • {OUTPUT_FASTA}")
print(f"  • {OUTPUT_CLUSTAL}")
print(f"  • {OUTPUT_STATS}")
print()
print("Next: Phase 3 - Conservation analysis and validation")
print("       Compare conserved positions against ground truth")
print()
