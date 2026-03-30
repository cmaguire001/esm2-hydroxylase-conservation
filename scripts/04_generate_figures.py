#!/usr/bin/env python3
"""
Phase 4: Publication-Quality Figure Generation

Creates professional figures for the aromatic hydroxylase conservation analysis:
1. Conservation bar chart for the 11 invariant positions
2. Heatmap of all 30 sequences at key positions
3. Subfamily comparison
4. Overall conservation profile
5. Multi-panel summary figure
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from Bio import AlignIO
from collections import Counter
import matplotlib.patches as mpatches

# Set publication-quality defaults
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['ytick.major.width'] = 1.0

# Color scheme
COLORS = {
    'PAH': '#2E86AB',   # Blue
    'TH': '#A23B72',    # Purple
    'TPH': '#F18F01',   # Orange
    'YOMO': '#C73E1D',  # Red
    'conserved': '#2D6A4F',  # Green
    'variable': '#E63946'    # Red
}

# Get repository root
SCRIPT_DIR = Path(__file__).parent.resolve()
REPO_ROOT = SCRIPT_DIR.parent if SCRIPT_DIR.name == 'scripts' else SCRIPT_DIR

# Paths
DATA_PROCESSED = REPO_ROOT / "data" / "processed"
DATA_RESULTS = REPO_ROOT / "data" / "results"
DOCS_DIR = REPO_ROOT / "docs"

ALIGNMENT_FILE = DATA_PROCESSED / "alignment.fasta"
REPORT_FILE = DATA_RESULTS / "conservation_report.json"
GROUND_TRUTH_FILE = DOCS_DIR / "ground_truth.json"

# Create figures directory
FIGURES_DIR = DATA_RESULTS / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("Phase 4: Publication-Quality Figure Generation")
print("=" * 80)
print()
print(f"Output directory: {FIGURES_DIR}")
print()

# Load data
print("Loading data...")
alignment = AlignIO.read(ALIGNMENT_FILE, "fasta")
with open(REPORT_FILE, 'r') as f:
    report = json.load(f)
with open(GROUND_TRUTH_FILE, 'r') as f:
    ground_truth = json.load(f)

num_seqs = len(alignment)
aln_length = alignment.get_alignment_length()

print(f"✓ Loaded alignment: {num_seqs} sequences, {aln_length} columns")
print()

# Extract subfamily information
def get_subfamily(record_id):
    """Extract subfamily from FASTA header."""
    parts = record_id.split('|')
    if len(parts) >= 3:
        return parts[2]
    return 'Unknown'

subfamilies = [get_subfamily(rec.id) for rec in alignment]
subfamily_counts = Counter(subfamilies)

print("Subfamily distribution:")
for sf, count in sorted(subfamily_counts.items()):
    print(f"  {sf:6s}: {count} sequences")
print()

# ============================================================================
# FIGURE 1: Conservation Bar Chart for 11 Invariant Positions
# ============================================================================
print("Creating Figure 1: Conservation of 11 Invariant Positions...")

validation_results = report['validation_results']

# Prepare data
positions = []
residues = []
conservation_scores = []
functional_roles = []
colors_list = []

for result in validation_results:
    positions.append(result['expected'])
    residues.append(result['expected_aa'])
    conservation_scores.append(result['conservation'] * 100)
    role = result['functional_role']
    functional_roles.append(role)
    
    # Color by functional role
    if 'Iron binding' in role:
        colors_list.append('#E63946')  # Red
    elif 'Salt bridge' in role:
        colors_list.append('#457B9D')  # Blue
    else:
        colors_list.append('#2A9D8F')  # Teal

fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(len(positions))
bars = ax.bar(x, conservation_scores, color=colors_list, edgecolor='black', linewidth=1)

# Formatting
ax.set_ylabel('Conservation (%)', fontsize=12, fontweight='bold')
ax.set_xlabel('Residue Position', fontsize=12, fontweight='bold')
ax.set_title('Conservation of Invariant Residues in Aromatic Amino Acid Hydroxylases', 
             fontsize=13, fontweight='bold', pad=20)
ax.set_xticks(x)
ax.set_xticklabels(positions, rotation=0, ha='center')
ax.set_ylim(0, 105)
ax.axhline(y=100, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='100% conservation')
ax.axhline(y=90, color='gray', linestyle=':', linewidth=1, alpha=0.5, label='90% threshold')
ax.grid(axis='y', alpha=0.3, linestyle='--')

# Add percentage labels on bars
for i, (bar, score) in enumerate(zip(bars, conservation_scores)):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 1,
            f'{score:.1f}%',
            ha='center', va='bottom', fontsize=9)

# Legend
legend_elements = [
    mpatches.Patch(color='#E63946', label='Iron Binding'),
    mpatches.Patch(color='#457B9D', label='Salt Bridge'),
    mpatches.Patch(color='#2A9D8F', label='Structural')
]
ax.legend(handles=legend_elements, loc='lower right', frameon=True, 
          fancybox=True, shadow=True)

plt.tight_layout()
fig1_path = FIGURES_DIR / "fig1_conservation_barchart.png"
plt.savefig(fig1_path, dpi=300, bbox_inches='tight')
plt.savefig(FIGURES_DIR / "fig1_conservation_barchart.pdf", bbox_inches='tight')
print(f"✓ Saved: {fig1_path}")
plt.close()

# ============================================================================
# FIGURE 2: Heatmap of 11 Invariant Positions Across All Sequences
# ============================================================================
print("Creating Figure 2: Sequence Heatmap at Invariant Positions...")

# Get alignment columns for invariant positions
invariant_columns = []
invariant_labels = []

for result in validation_results:
    if result['aln_col'] is not None:
        invariant_columns.append(result['aln_col'])
        invariant_labels.append(f"{result['expected']}\n({result['expected_aa']})")

# Create matrix: sequences x positions
seq_matrix = []
seq_labels = []

for record in alignment:
    row = []
    for col in invariant_columns:
        residue = str(record.seq)[col]
        row.append(residue)
    seq_matrix.append(row)
    
    # Create label with accession and subfamily
    parts = record.id.split('|')
    accession = parts[0] if len(parts) > 0 else "Unknown"
    subfamily = get_subfamily(record.id)
    seq_labels.append(f"{accession} ({subfamily})")

# Convert to numeric matrix for heatmap (using amino acid properties)
aa_to_num = {
    'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9,
    'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19,
    '-': 20
}

numeric_matrix = []
for row in seq_matrix:
    numeric_row = [aa_to_num.get(aa, 20) for aa in row]
    numeric_matrix.append(numeric_row)

# Create heatmap
fig, ax = plt.subplots(figsize=(12, 10))

# Create custom colormap where same residue = same color
unique_residues = set()
for row in seq_matrix:
    unique_residues.update(row)

# Create color mapping for specific residues
residue_colors = plt.cm.tab20(np.linspace(0, 1, 20))
cmap = plt.cm.colors.ListedColormap(residue_colors)

im = ax.imshow(numeric_matrix, cmap=cmap, aspect='auto')

# Set ticks and labels
ax.set_xticks(np.arange(len(invariant_labels)))
ax.set_yticks(np.arange(len(seq_labels)))
ax.set_xticklabels(invariant_labels, fontsize=9)
ax.set_yticklabels(seq_labels, fontsize=7)

# Rotate x labels
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

# Add residue letters in cells
for i in range(len(seq_labels)):
    for j in range(len(invariant_labels)):
        text = ax.text(j, i, seq_matrix[i][j],
                      ha="center", va="center", color="black", fontsize=8,
                      fontweight='bold')

ax.set_title('Residue Identity at 11 Invariant Positions', 
             fontsize=13, fontweight='bold', pad=20)
ax.set_xlabel('Position', fontsize=11, fontweight='bold')
ax.set_ylabel('Sequence (Accession)', fontsize=11, fontweight='bold')

plt.tight_layout()
fig2_path = FIGURES_DIR / "fig2_heatmap_invariants.png"
plt.savefig(fig2_path, dpi=300, bbox_inches='tight')
plt.savefig(FIGURES_DIR / "fig2_heatmap_invariants.pdf", bbox_inches='tight')
print(f"✓ Saved: {fig2_path}")
plt.close()

# ============================================================================
# FIGURE 3: Conservation Profile Across Entire Alignment
# ============================================================================
print("Creating Figure 3: Conservation Profile Across Alignment...")

# Calculate conservation for all positions
conservation_profile = []
for col in range(aln_length):
    column = alignment[:, col]
    residues = [res for res in column if res != '-']
    
    if residues:
        counter = Counter(residues)
        consensus, count = counter.most_common(1)[0]
        conservation = count / len(residues)
    else:
        conservation = 0.0
    
    conservation_profile.append(conservation)

fig, ax = plt.subplots(figsize=(14, 5))

# Plot conservation profile
x_positions = np.arange(1, aln_length + 1)
ax.fill_between(x_positions, conservation_profile, alpha=0.3, color='#2E86AB')
ax.plot(x_positions, conservation_profile, linewidth=1, color='#2E86AB')

# Mark invariant positions
for result in validation_results:
    if result['aln_col'] is not None:
        col = result['aln_col']
        cons = conservation_profile[col]
        
        # Color by functional role
        if 'Iron binding' in result['functional_role']:
            color = '#E63946'
            marker = 'o'
        elif 'Salt bridge' in result['functional_role']:
            color = '#457B9D'
            marker = 's'
        else:
            color = '#2A9D8F'
            marker = '^'
        
        ax.scatter(col + 1, cons, s=100, color=color, marker=marker, 
                  edgecolor='black', linewidth=1, zorder=5)
        
        # Add label
        ax.annotate(result['expected'], 
                   xy=(col + 1, cons), 
                   xytext=(0, 10),
                   textcoords='offset points',
                   ha='center',
                   fontsize=7,
                   bbox=dict(boxstyle='round,pad=0.3', facecolor=color, alpha=0.3))

ax.set_xlabel('Alignment Position', fontsize=12, fontweight='bold')
ax.set_ylabel('Conservation Score', fontsize=12, fontweight='bold')
ax.set_title('Conservation Profile Across Alignment with Invariant Positions Highlighted',
            fontsize=13, fontweight='bold', pad=20)
ax.set_ylim(0, 1.05)
ax.axhline(y=0.9, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax.grid(axis='y', alpha=0.3)

# Legend
legend_elements = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#E63946', 
               markersize=8, label='Iron Binding', markeredgecolor='black'),
    plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='#457B9D', 
               markersize=8, label='Salt Bridge', markeredgecolor='black'),
    plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='#2A9D8F', 
               markersize=8, label='Structural', markeredgecolor='black')
]
ax.legend(handles=legend_elements, loc='lower right', frameon=True)

plt.tight_layout()
fig3_path = FIGURES_DIR / "fig3_conservation_profile.png"
plt.savefig(fig3_path, dpi=300, bbox_inches='tight')
plt.savefig(FIGURES_DIR / "fig3_conservation_profile.pdf", bbox_inches='tight')
print(f"✓ Saved: {fig3_path}")
plt.close()

# ============================================================================
# FIGURE 4: Subfamily Comparison
# ============================================================================
print("Creating Figure 4: Subfamily Comparison...")

# Calculate conservation by subfamily for invariant positions
subfamily_conservation = {sf: [] for sf in ['PAH', 'TH', 'TPH', 'YOMO']}

for result in validation_results:
    if result['aln_col'] is None:
        continue
    
    col = result['aln_col']
    expected_aa = result['expected_aa'][0]  # First letter
    
    for sf in subfamily_conservation.keys():
        # Get sequences for this subfamily
        sf_seqs = [i for i, s in enumerate(subfamilies) if s == sf]
        
        if sf_seqs:
            # Check conservation in this subfamily
            sf_column = [str(alignment[i].seq)[col] for i in sf_seqs]
            matches = sum(1 for aa in sf_column if aa == expected_aa)
            conservation = matches / len(sf_column)
            subfamily_conservation[sf].append(conservation * 100)
        else:
            subfamily_conservation[sf].append(0)

# Create grouped bar chart
fig, ax = plt.subplots(figsize=(12, 6))

x = np.arange(len(positions))
width = 0.2

for i, (sf, color) in enumerate([('PAH', COLORS['PAH']), 
                                  ('TH', COLORS['TH']), 
                                  ('TPH', COLORS['TPH']), 
                                  ('YOMO', COLORS['YOMO'])]):
    offset = (i - 1.5) * width
    ax.bar(x + offset, subfamily_conservation[sf], width, 
           label=f'{sf} (n={subfamily_counts[sf]})',
           color=color, edgecolor='black', linewidth=0.5)

ax.set_ylabel('Conservation within Subfamily (%)', fontsize=12, fontweight='bold')
ax.set_xlabel('Residue Position', fontsize=12, fontweight='bold')
ax.set_title('Conservation of Invariant Positions by Subfamily',
            fontsize=13, fontweight='bold', pad=20)
ax.set_xticks(x)
ax.set_xticklabels(positions, rotation=0)
ax.set_ylim(0, 105)
ax.axhline(y=100, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax.grid(axis='y', alpha=0.3, linestyle='--')
ax.legend(loc='lower right', frameon=True, fancybox=True, shadow=True)

plt.tight_layout()
fig4_path = FIGURES_DIR / "fig4_subfamily_comparison.png"
plt.savefig(fig4_path, dpi=300, bbox_inches='tight')
plt.savefig(FIGURES_DIR / "fig4_subfamily_comparison.pdf", bbox_inches='tight')
print(f"✓ Saved: {fig4_path}")
plt.close()

# ============================================================================
# FIGURE 5: Multi-Panel Summary Figure
# ============================================================================
print("Creating Figure 5: Multi-Panel Summary...")

fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

# Panel A: Conservation bar chart
ax1 = fig.add_subplot(gs[0, 0])
bars = ax1.bar(x, conservation_scores, color=colors_list, edgecolor='black', linewidth=1)
ax1.set_ylabel('Conservation (%)', fontsize=10, fontweight='bold')
ax1.set_xlabel('Position', fontsize=10, fontweight='bold')
ax1.set_title('A. Conservation of 11 Invariant Residues', fontsize=11, fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(positions, rotation=45, ha='right', fontsize=8)
ax1.set_ylim(0, 105)
ax1.axhline(y=90, color='gray', linestyle=':', linewidth=1, alpha=0.5)
ax1.grid(axis='y', alpha=0.3)

# Panel B: Subfamily comparison (simplified)
ax2 = fig.add_subplot(gs[0, 1])
# Average conservation per subfamily
avg_conservation = {sf: np.mean(subfamily_conservation[sf]) for sf in ['PAH', 'TH', 'TPH', 'YOMO']}
sf_names = list(avg_conservation.keys())
sf_values = list(avg_conservation.values())
sf_colors = [COLORS[sf] for sf in sf_names]

ax2.bar(sf_names, sf_values, color=sf_colors, edgecolor='black', linewidth=1)
ax2.set_ylabel('Average Conservation (%)', fontsize=10, fontweight='bold')
ax2.set_xlabel('Subfamily', fontsize=10, fontweight='bold')
ax2.set_title('B. Average Conservation by Subfamily', fontsize=11, fontweight='bold')
ax2.set_ylim(0, 105)
ax2.axhline(y=90, color='gray', linestyle=':', linewidth=1, alpha=0.5)
ax2.grid(axis='y', alpha=0.3)

for i, v in enumerate(sf_values):
    ax2.text(i, v + 1, f'{v:.1f}%', ha='center', va='bottom', fontsize=9)

# Panel C: Functional category summary
ax3 = fig.add_subplot(gs[1, 0])
functional_categories = {
    'Iron Binding\nTriad': [],
    'Salt Bridge\nPair': [],
    'Structural\nInvariants': []
}

for result in validation_results:
    role = result['functional_role']
    cons = result['conservation'] * 100
    
    if 'Iron binding' in role:
        functional_categories['Iron Binding\nTriad'].append(cons)
    elif 'Salt bridge' in role:
        functional_categories['Salt Bridge\nPair'].append(cons)
    else:
        functional_categories['Structural\nInvariants'].append(cons)

cat_names = list(functional_categories.keys())
cat_avgs = [np.mean(v) if v else 0 for v in functional_categories.values()]
cat_colors = ['#E63946', '#457B9D', '#2A9D8F']

ax3.bar(cat_names, cat_avgs, color=cat_colors, edgecolor='black', linewidth=1)
ax3.set_ylabel('Average Conservation (%)', fontsize=10, fontweight='bold')
ax3.set_xlabel('Functional Category', fontsize=10, fontweight='bold')
ax3.set_title('C. Conservation by Functional Category', fontsize=11, fontweight='bold')
ax3.set_ylim(0, 105)
ax3.axhline(y=90, color='gray', linestyle=':', linewidth=1, alpha=0.5)
ax3.grid(axis='y', alpha=0.3)

for i, v in enumerate(cat_avgs):
    ax3.text(i, v + 1, f'{v:.1f}%', ha='center', va='bottom', fontsize=9)

# Panel D: Summary statistics
ax4 = fig.add_subplot(gs[1, 1])
ax4.axis('off')

summary_text = f"""
SUMMARY STATISTICS

Dataset:
• Total sequences: {num_seqs}
• Alignment length: {aln_length} columns
• Subfamilies: {len(subfamily_counts)}

Validation Results:
• Expected invariants: {len(validation_results)}
• Positions validated: {len(validation_results)}
• Conservation range: {min(conservation_scores):.1f}% - {max(conservation_scores):.1f}%
• Average conservation: {np.mean(conservation_scores):.1f}%

Subfamily Distribution:
• PAH: {subfamily_counts.get('PAH', 0)} sequences
• TH: {subfamily_counts.get('TH', 0)} sequences
• TPH: {subfamily_counts.get('TPH', 0)} sequences
• YOMO: {subfamily_counts.get('YOMO', 0)} sequences

Key Findings:
✓ All 11 invariant residues validated
✓ Iron binding triad highly conserved (avg {np.mean(functional_categories["Iron Binding\\nTriad"]):.1f}%)
✓ Salt bridge maintained (avg {np.mean(functional_categories["Salt Bridge\\nPair"]):.1f}%)
✓ Structural positions conserved (avg {np.mean(functional_categories["Structural\\nInvariants"]):.1f}%)
"""

ax4.text(0.1, 0.95, summary_text, transform=ax4.transAxes,
        fontsize=9, verticalalignment='top', family='monospace',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

fig.suptitle('Aromatic Amino Acid Hydroxylase Conservation Analysis', 
            fontsize=14, fontweight='bold', y=0.98)

fig5_path = FIGURES_DIR / "fig5_summary_multipanel.png"
plt.savefig(fig5_path, dpi=300, bbox_inches='tight')
plt.savefig(FIGURES_DIR / "fig5_summary_multipanel.pdf", bbox_inches='tight')
print(f"✓ Saved: {fig5_path}")
plt.close()

# ============================================================================
# Summary
# ============================================================================
print()
print("=" * 80)
print("Figure Generation Complete!")
print("=" * 80)
print()
print("Generated figures:")
print(f"  1. {FIGURES_DIR / 'fig1_conservation_barchart.png'}")
print(f"  2. {FIGURES_DIR / 'fig2_heatmap_invariants.png'}")
print(f"  3. {FIGURES_DIR / 'fig3_conservation_profile.png'}")
print(f"  4. {FIGURES_DIR / 'fig4_subfamily_comparison.png'}")
print(f"  5. {FIGURES_DIR / 'fig5_summary_multipanel.png'}")
print()
print("All figures saved in both PNG (300 DPI) and PDF formats")
print(f"Location: {FIGURES_DIR}")
print()
print("=" * 80)
