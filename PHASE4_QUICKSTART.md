# Phase 4: Publication Figure Generation - Quick Start

## What Phase 4 Creates

Five publication-ready figures (PNG at 300 DPI + PDF):

1. **Conservation Bar Chart** - Shows conservation % for all 11 invariant positions, color-coded by function
2. **Heatmap** - All 30 sequences × 11 positions, showing exact residues at each position
3. **Conservation Profile** - Full alignment conservation with invariant positions highlighted
4. **Subfamily Comparison** - Conservation broken down by PAH/TH/TPH/YOMO
5. **Multi-Panel Summary** - Combined figure with all key results

## Prerequisites

Matplotlib and seaborn should already be installed from requirements.txt, but verify:

```bash
pip install matplotlib seaborn numpy
```

## Run Phase 4

```bash
cd /workspaces/esm2-hydroxylase-conservation

# Copy script to scripts directory
mv 04_generate_figures.py scripts/

# Run figure generation
cd scripts
python 04_generate_figures.py
```

**Expected runtime:** 10-20 seconds

## Output

All figures saved to: `data/results/figures/`

```
data/results/figures/
├── fig1_conservation_barchart.png (+ .pdf)
├── fig2_heatmap_invariants.png (+ .pdf)
├── fig3_conservation_profile.png (+ .pdf)
├── fig4_subfamily_comparison.png (+ .pdf)
└── fig5_summary_multipanel.png (+ .pdf)
```

## Figure Descriptions

### Figure 1: Conservation Bar Chart
- **Use for:** Showing overall validation success
- **Shows:** Conservation % for each of 11 positions
- **Colors:** Red = Iron binding, Blue = Salt bridge, Teal = Structural
- **Key finding:** All positions >89% conserved

### Figure 2: Heatmap at Invariant Positions
- **Use for:** Showing sequence-level detail
- **Shows:** Actual residues for all 30 sequences at 11 positions
- **Reveals:** Which sequences have variants, subfamily patterns

### Figure 3: Conservation Profile
- **Use for:** Genome-wide view of conservation
- **Shows:** Conservation across all 809 alignment columns
- **Highlights:** Your 11 invariant positions stand out as peaks
- **Markers:** Circle = Iron binding, Square = Salt bridge, Triangle = Structural

### Figure 4: Subfamily Comparison
- **Use for:** Comparing PAH vs TH vs TPH vs YOMO
- **Shows:** Grouped bars for each subfamily at each position
- **Reveals:** Which subfamilies are more/less conserved

### Figure 5: Multi-Panel Summary (Publication Main Figure)
- **Use for:** Manuscript main figure or poster
- **Contains:**
  - Panel A: Conservation bar chart
  - Panel B: Average conservation by subfamily
  - Panel C: Conservation by functional category
  - Panel D: Summary statistics
- **Size:** 16×10 inches, perfect for full-page figure

## Customization

To adjust figure appearance, edit these parameters in the script:

```python
# DPI (resolution)
plt.rcParams['figure.dpi'] = 300  # Change to 600 for ultra-high-res

# Font sizes
plt.rcParams['font.size'] = 10  # Adjust for readability

# Colors
COLORS = {
    'PAH': '#2E86AB',   # Change hex codes for different colors
    'TH': '#A23B72',
    ...
}
```

## Viewing Figures

In Codespaces:
1. Navigate to `data/results/figures/` in file explorer
2. Click any `.png` file to preview
3. Download for use in presentations/papers

Or use command line:
```bash
# List all figures
ls -lh data/results/figures/

# Open in viewer (if X11 forwarding enabled)
eog data/results/figures/fig5_summary_multipanel.png
```

## Publication Tips

**For manuscripts:**
- Use PDF versions (vector graphics, scalable)
- Figure 5 is perfect for main text
- Figures 1-4 as supplementary

**For presentations:**
- Use PNG at 300 DPI
- Figure 1 or 5 for overview slides
- Figure 2 or 4 for detailed discussion

**For posters:**
- Consider increasing DPI to 600
- Figure 5 works great as centerpiece
- High contrast colors print well

## Troubleshooting

**"No module named 'matplotlib'"**
```bash
pip install matplotlib seaborn
```

**Fonts look wrong**
- Script defaults to Arial/DejaVu Sans
- Should work on most systems
- Edit font.sans-serif if needed

**Figures too small/large**
- Edit figsize parameters in script
- Example: `figsize=(10, 6)` → `figsize=(12, 8)`

**Colors not showing**
- Check if running in headless mode
- Try different backend: `plt.switch_backend('Agg')`

---

**Ready to generate?** Run the script and you'll have publication-ready figures in ~20 seconds! 📊
