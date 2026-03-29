#!/usr/bin/env python3
"""
Phase 3: Conservation Analysis and Validation

This script:
1. Maps ground truth positions to alignment columns
2. Identifies the catalytic domain region
3. Calculates conservation in the functional core
4. Validates against the 11 invariant residues from 2011

The challenge: Position numbers in ground truth (e.g., H-285) refer to 
positions in a reference sequence (likely human PAH catalytic domain),
not alignment column numbers.
"""

import json
from pathlib import Path
from collections import Counter
from Bio import AlignIO, SeqIO
import re

# Get repository root
SCRIPT_DIR = Path(__file__).parent.resolve()
REPO_ROOT = SCRIPT_DIR.parent if SCRIPT_DIR.name == 'scripts' else SCRIPT_DIR

# Paths
DATA_PROCESSED = REPO_ROOT / "data" / "processed"
DATA_RESULTS = REPO_ROOT / "data" / "results"
DATA_RESULTS.mkdir(parents=True, exist_ok=True)

DOCS_DIR = REPO_ROOT / "docs"

ALIGNMENT_FILE = DATA_PROCESSED / "alignment.fasta"
STATS_FILE = DATA_PROCESSED / "alignment_stats.json"
GROUND_TRUTH_FILE = DOCS_DIR / "ground_truth.json"

OUTPUT_REPORT = DATA_RESULTS / "conservation_report.json"
OUTPUT_VALIDATION = DATA_RESULTS / "validation_report.txt"

print("=" * 80)
print("Phase 3: Conservation Analysis & Validation")
print("=" * 80)
print()

# Load alignment
if not ALIGNMENT_FILE.exists():
    print(f"❌ ERROR: Alignment file not found: {ALIGNMENT_FILE}")
    print("Run Phase 2 (02_align_sequences.py) first!")
    exit(1)

alignment = AlignIO.read(ALIGNMENT_FILE, "fasta")
num_seqs = len(alignment)
aln_length = alignment.get_alignment_length()

print(f"Loaded alignment: {num_seqs} sequences, {aln_length} columns")
print()

# Load ground truth
if not GROUND_TRUTH_FILE.exists():
    print(f"⚠️  WARNING: Ground truth file not found: {GROUND_TRUTH_FILE}")
    print("Proceeding with conservation analysis only...")
    ground_truth = None
else:
    with open(GROUND_TRUTH_FILE, 'r') as f:
        ground_truth = json.load(f)
    
    invariants = ground_truth['invariant_residues']['residues']
    print(f"Loaded ground truth: {len(invariants)} invariant residues expected")
    print()

# Find reference sequence (human PAH)
print("=" * 80)
print("Finding Reference Sequence")
print("=" * 80)
print()

reference_record = None
reference_index = None

for i, record in enumerate(alignment):
    # Look for human PAH (P00439 or contains "HUMAN" and "PAH")
    if 'P00439' in record.id or ('HUMAN' in record.id.upper() and 'PAH' in record.id.upper()):
        reference_record = record
        reference_index = i
        print(f"✓ Found reference: {record.id}")
        print(f"  Position in alignment: {i}")
        break

if not reference_record:
    # Fallback: use first PAH sequence
    for i, record in enumerate(alignment):
        if 'PAH' in record.id.upper():
            reference_record = record
            reference_index = i
            print(f"⚠️  Using first PAH sequence as reference: {record.id}")
            break

if not reference_record:
    print("❌ ERROR: No suitable reference sequence found!")
    print("Need a PAH sequence to map ground truth positions.")
    exit(1)

print()

# Map reference sequence positions to alignment columns
print("=" * 80)
print("Mapping Reference Positions to Alignment Columns")
print("=" * 80)
print()

ref_seq = str(reference_record.seq)
ref_seq_nogaps = ref_seq.replace('-', '')

print(f"Reference sequence length (no gaps): {len(ref_seq_nogaps)} aa")
print(f"Alignment length: {len(ref_seq)} columns")
print()

# Create mapping: reference position -> alignment column
ref_pos_to_aln_col = {}
aln_col_to_ref_pos = {}

ref_pos = 0  # Position in ungapped reference (1-indexed after increment)
for aln_col in range(len(ref_seq)):
    if ref_seq[aln_col] != '-':
        ref_pos += 1
        ref_pos_to_aln_col[ref_pos] = aln_col
        aln_col_to_ref_pos[aln_col] = ref_pos

print(f"Created position mapping for {len(ref_pos_to_aln_col)} reference positions")
print()

# Example mappings
print("Example position mappings:")
example_positions = [1, 50, 100, 200, 285, 290, 330]
for pos in example_positions:
    if pos in ref_pos_to_aln_col:
        aln_col = ref_pos_to_aln_col[pos]
        residue = ref_seq[aln_col]
        print(f"  Ref pos {pos:3d} -> Alignment col {aln_col:4d} ({residue})")
print()

# Calculate conservation for all positions
print("=" * 80)
print("Calculating Conservation Scores")
print("=" * 80)
print()

conservation_scores = []
highly_conserved = []  # >90% conservation
fully_conserved = []   # 100% conservation

for aln_col in range(aln_length):
    column = alignment[:, aln_col]
    
    # Get residues (excluding gaps)
    residues = [res for res in column if res != '-']
    
    if not residues:
        conservation = 0.0
        consensus = '-'
        identity = 0
    else:
        counter = Counter(residues)
        consensus, count = counter.most_common(1)[0]
        
        # Conservation = fraction of sequences with consensus (among non-gap)
        conservation = count / len(residues)
        
        # Identity = fraction of ALL sequences with consensus
        identity = count / num_seqs
        
        if conservation == 1.0 and len(residues) == num_seqs:
            fully_conserved.append({
                'aln_col': aln_col,
                'ref_pos': aln_col_to_ref_pos.get(aln_col, None),
                'residue': consensus,
                'conservation': conservation
            })
        elif conservation >= 0.9:
            highly_conserved.append({
                'aln_col': aln_col,
                'ref_pos': aln_col_to_ref_pos.get(aln_col, None),
                'residue': consensus,
                'conservation': conservation
            })
    
    conservation_scores.append({
        'aln_col': aln_col,
        'ref_pos': aln_col_to_ref_pos.get(aln_col, None),
        'consensus': consensus,
        'conservation': conservation,
        'num_sequences': len(residues),
        'gaps': num_seqs - len(residues)
    })

print(f"✓ Fully conserved (100%): {len(fully_conserved)} positions")
print(f"✓ Highly conserved (≥90%): {len(highly_conserved)} positions")
print()

# Display conserved positions
if fully_conserved:
    print("Fully conserved positions:")
    for pos in fully_conserved:
        ref_info = f"Ref:{pos['ref_pos']:3d}" if pos['ref_pos'] else "Ref:N/A"
        print(f"  Aln col {pos['aln_col']:4d} ({ref_info}): {pos['residue']}")
print()

# Validate against ground truth
if ground_truth:
    print("=" * 80)
    print("Validating Against Ground Truth")
    print("=" * 80)
    print()
    
    validation_results = []
    
    for invariant in invariants:
        expected_aa = invariant['amino_acid']
        expected_pos = invariant['position']
        expected_code = invariant['code']
        functional_role = invariant.get('functional_role', 'Unknown')
        
        # Map to alignment column
        if expected_pos in ref_pos_to_aln_col:
            aln_col = ref_pos_to_aln_col[expected_pos]
            
            # Get conservation at this position
            col_stats = conservation_scores[aln_col]
            consensus = col_stats['consensus']
            conservation = col_stats['conservation']
            
            # Check if it matches
            matches = (consensus == expected_aa[0])  # First letter of amino acid name
            
            # Get actual column
            column = alignment[:, aln_col]
            residues = [res for res in column if res != '-']
            unique_residues = set(residues)
            
            validation_results.append({
                'expected': expected_code,
                'expected_aa': expected_aa,
                'expected_pos': expected_pos,
                'aln_col': aln_col,
                'observed_consensus': consensus,
                'observed_residues': list(unique_residues),
                'conservation': conservation,
                'matches': matches,
                'functional_role': functional_role
            })
            
            status = "✓" if matches and conservation == 1.0 else "✗" if not matches else "~"
            print(f"{status} {expected_code:8s} (Ref pos {expected_pos:3d} -> Aln col {aln_col:4d})")
            print(f"    Expected: {expected_aa:3s}  Observed: {consensus} ({conservation*100:.1f}% conserved)")
            print(f"    Role: {functional_role}")
            if not matches:
                print(f"    ⚠️  MISMATCH! Found: {unique_residues}")
            print()
        else:
            print(f"✗ {expected_code:8s} (Ref pos {expected_pos:3d})")
            print(f"    ❌ Position not found in reference sequence!")
            print(f"    (Reference sequence may be truncated)")
            print()
            
            validation_results.append({
                'expected': expected_code,
                'expected_aa': expected_aa,
                'expected_pos': expected_pos,
                'aln_col': None,
                'observed_consensus': None,
                'conservation': 0.0,
                'matches': False,
                'functional_role': functional_role,
                'error': 'Position not in reference sequence'
            })
    
    # Summary
    validated = sum(1 for r in validation_results if r['matches'] and r['conservation'] == 1.0)
    partial = sum(1 for r in validation_results if r['matches'] and r['conservation'] < 1.0)
    failed = sum(1 for r in validation_results if not r['matches'])
    missing = sum(1 for r in validation_results if r['aln_col'] is None)
    
    print("=" * 80)
    print("VALIDATION SUMMARY")
    print("=" * 80)
    print(f"Expected invariant residues: {len(invariants)}")
    print(f"Fully validated (100% conserved): {validated}")
    print(f"Partially validated (present but <100%): {partial}")
    print(f"Failed (different residue): {failed}")
    print(f"Missing (position not in reference): {missing}")
    print()
    
    # Save validation report
    report = {
        'alignment_stats': {
            'num_sequences': num_seqs,
            'alignment_length': aln_length,
            'fully_conserved_positions': len(fully_conserved),
            'highly_conserved_positions': len(highly_conserved)
        },
        'ground_truth': {
            'expected_invariants': len(invariants),
            'validated': validated,
            'partial': partial,
            'failed': failed,
            'missing': missing
        },
        'validation_results': validation_results,
        'fully_conserved_positions': fully_conserved,
        'highly_conserved_positions': highly_conserved
    }
    
    with open(OUTPUT_REPORT, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"✓ Detailed report saved: {OUTPUT_REPORT}")
    
    # Generate text summary
    with open(OUTPUT_VALIDATION, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("AROMATIC AMINO ACID HYDROXYLASE CONSERVATION ANALYSIS\n")
        f.write("Validation Report\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Analysis Date: 2026-03-29\n")
        f.write(f"Sequences Analyzed: {num_seqs}\n")
        f.write(f"Alignment Length: {aln_length} columns\n\n")
        
        f.write("GROUND TRUTH VALIDATION\n")
        f.write("-" * 80 + "\n\n")
        
        f.write(f"Expected Invariant Residues: {len(invariants)}\n")
        f.write(f"Fully Validated: {validated}\n")
        f.write(f"Partially Validated: {partial}\n")
        f.write(f"Failed: {failed}\n")
        f.write(f"Missing: {missing}\n\n")
        
        f.write("DETAILED RESULTS\n")
        f.write("-" * 80 + "\n\n")
        
        for result in validation_results:
            status = "PASS" if result['matches'] and result['conservation'] == 1.0 else "PARTIAL" if result['matches'] else "FAIL"
            f.write(f"{result['expected']:10s} [{status:8s}]\n")
            f.write(f"  Expected: {result['expected_aa']} at position {result['expected_pos']}\n")
            if result['aln_col'] is not None:
                f.write(f"  Observed: {result['observed_consensus']} ({result['conservation']*100:.1f}% conserved)\n")
                f.write(f"  Alignment column: {result['aln_col']}\n")
            else:
                f.write(f"  ERROR: Position not found in reference sequence\n")
            f.write(f"  Function: {result['functional_role']}\n\n")
    
    print(f"✓ Text summary saved: {OUTPUT_VALIDATION}")

else:
    # No ground truth - just report conservation
    report = {
        'alignment_stats': {
            'num_sequences': num_seqs,
            'alignment_length': aln_length,
            'fully_conserved_positions': len(fully_conserved),
            'highly_conserved_positions': len(highly_conserved)
        },
        'fully_conserved_positions': fully_conserved,
        'highly_conserved_positions': highly_conserved,
        'all_positions': conservation_scores
    }
    
    with open(OUTPUT_REPORT, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"✓ Conservation report saved: {OUTPUT_REPORT}")

print()
print("=" * 80)
print("Phase 3 Complete!")
print("=" * 80)
print()
