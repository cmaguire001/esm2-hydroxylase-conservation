# Understanding Phase 2 Results & Phase 3 Strategy

## What Happened in Phase 2?

We aligned 30 sequences and found only **5 fully conserved positions**, but expected **11 invariant residues** from the 2011 analysis.

### Why the Discrepancy?

**The Core Issue: Domain Architecture**

Your sequences include different regions:
- **Full-length proteins** (with regulatory domains, targeting signals)
- **Different isoforms** (TH has long N-terminal regulatory regions)
- **Variable termini** (different start/end points)

Evidence:
- Alignment: 809 columns
- Average sequence: 468 amino acids
- Gap percentage: 42% (very high)
- Alignment preview shows dashes at the start

### What This Means

The **position numbers in your ground truth** (H-285, H-290, E-330, etc.) refer to positions in the **catalytic domain** of a reference sequence (likely human PAH), NOT to alignment column numbers.

**Example:**
- Ground truth says "H-285" (Histidine at position 285)
- This is position 285 in **human PAH catalytic domain**
- In the full alignment, this might be column 412 (after N-terminal gaps)
- Different sequences have different N-terminal lengths!

### The 2011 Analysis Probably Used:

1. **Only catalytic domains** (trimmed to conserved core)
2. **Consistent domain boundaries** (all sequences start/end at same point)
3. **Position numbering from one reference** (likely human PAH)

## Phase 3 Solution: Position Mapping

The Phase 3 script solves this by:

### Step 1: Find Reference Sequence
- Identifies human PAH (P00439) in the alignment
- Uses it as the numbering reference

### Step 2: Create Position Map
- Maps reference sequence positions → alignment columns
- Example: "Position 285 in PAH" → "Column 412 in alignment"

### Step 3: Look Up Ground Truth
For each of your 11 invariants:
1. Take the position number (e.g., 285)
2. Find the corresponding alignment column (e.g., 412)
3. Check conservation at that column
4. Validate the amino acid matches

### Step 4: Report
- How many of the 11 invariants are conserved?
- Which ones match expectations?
- Which positions have variants?

## Expected Outcomes

### Best Case Scenario
All 11 positions show 100% conservation when we look at the right alignment columns:
- ✓ H-285, H-290, E-330 (iron binding triad)
- ✓ R-270, D-282 (salt bridge)
- ✓ S-349, P-281, G-289, G-344, P-225, E-353 (structural)

### Realistic Scenario
Most positions conserved, some might show:
- **High conservation (90-95%)** - mostly conserved with 1-2 variants
- **Conservative substitutions** - similar amino acids (e.g., Asp→Glu)
- **Subfamily differences** - YOMO might differ from PAH/TH/TPH

### What If Some Still Don't Match?

Possible reasons:
1. **Different reference sequence** - 2011 used different numbering
2. **Sequence quality** - Some UniProt entries might have errors
3. **Isoforms** - Different splice variants
4. **Subfamily divergence** - YOMO is more distant evolutionarily

## Running Phase 3

```bash
# Upload 03_analyze_conservation.py to your Codespace
# Move to scripts directory
mv 03_analyze_conservation.py scripts/

# Run it
cd scripts
python 03_analyze_conservation.py
```

**Expected output:**
```
Finding Reference Sequence
✓ Found reference: P00439|PH4H_HUMAN|...

Mapping Reference Positions to Alignment Columns
Created position mapping for 452 reference positions

Validating Against Ground Truth
✓ H-285 (Ref pos 285 -> Aln col 412)
  Expected: His  Observed: H (100.0% conserved)
  Role: Iron binding
...
```

## Success Metrics

**Full validation:**
- All 11 positions found in reference
- All 11 show expected amino acid
- High conservation (≥90%) for all

**Partial validation:**
- 8-10 positions fully conserved
- 1-3 positions show variants or conservative substitutions
- Still biologically meaningful!

## What We Learn

This is a **classic bioinformatics challenge**:
- Sequence databases contain diverse data (full-length, domains, isoforms)
- Position numbering is reference-dependent
- Domain boundaries matter for conservation analysis
- Need to map between coordinate systems

Your 2011 analysis was probably more manually curated (trimming sequences to just the catalytic core), while our automated approach uses full UniProt entries.

Both approaches are valid - manual curation gives cleaner results, automation is faster and reproducible!

---

**Ready for Phase 3?** Upload the script and let's see how many of your 11 invariants we can validate! 🧬
