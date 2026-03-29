#!/usr/bin/env python3
"""
Fetch aromatic amino acid hydroxylase sequences from UniProt.
Recreates 2011 analysis with representative sequences from 4 subfamilies:
- Phenylalanine hydroxylase (PAH)
- Tyrosine hydroxylase (TH)
- Tryptophan hydroxylase (TPH)
- Tyrosine monooxygenase

Author: Recreating 2011 poster analysis
Date: 2026
"""

import requests
import time
from pathlib import Path
from typing import List, Dict
import json

# Create output directory
DATA_DIR = Path("data/raw")
DATA_DIR.mkdir(parents=True, exist_ok=True)

# UniProt REST API base URL
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"

# Define search queries for each subfamily
SUBFAMILIES = {
    "PAH": {
        "name": "Phenylalanine hydroxylase",
        "query": "gene:PAH AND reviewed:true",
        "max_seqs": 8
    },
    "TH": {
        "name": "Tyrosine hydroxylase", 
        "query": "gene:TH AND reviewed:true",
        "max_seqs": 8
    },
    "TPH": {
        "name": "Tryptophan hydroxylase",
        "query": "(gene:TPH OR gene:TPH1 OR gene:TPH2) AND reviewed:true",
        "max_seqs": 8
    },
    "YOMO": {
        "name": "Tyrosine monooxygenase",
        # YOMO is often found in bacteria - using EC number for aromatic amino acid hydroxylases
        "query": "ec:1.14.16.* AND reviewed:true AND (tyrosine OR phenylalanine OR tryptophan) AND hydroxylase",
        "max_seqs": 6
    }
}

def fetch_sequences(query: str, max_results: int = 10) -> List[Dict]:
    """
    Fetch sequences from UniProt using REST API.
    
    Args:
        query: UniProt search query
        max_results: Maximum number of sequences to retrieve
        
    Returns:
        List of sequence records with metadata
    """
    params = {
        "query": query,
        "format": "json",
        "size": max_results,
        "fields": "accession,id,gene_names,organism_name,sequence,protein_name,length"
    }
    
    print(f"  Querying: {query}")
    print(f"  URL: {UNIPROT_API}")
    
    try:
        response = requests.get(UNIPROT_API, params=params, timeout=30)
        response.raise_for_status()
        
        data = response.json()
        results = data.get("results", [])
        
        print(f"  Found {len(results)} sequences")
        return results
        
    except requests.exceptions.RequestException as e:
        print(f"  ERROR: {e}")
        return []

def parse_sequence_record(record: Dict) -> Dict:
    """
    Parse UniProt JSON record into simplified format.
    
    Args:
        record: Raw UniProt JSON record
        
    Returns:
        Simplified sequence dictionary
    """
    return {
        "accession": record.get("primaryAccession", ""),
        "id": record.get("uniProtkbId", ""),
        "gene": record.get("genes", [{}])[0].get("geneName", {}).get("value", "") if record.get("genes") else "",
        "organism": record.get("organism", {}).get("scientificName", ""),
        "protein_name": record.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""),
        "length": record.get("sequence", {}).get("length", 0),
        "sequence": record.get("sequence", {}).get("value", "")
    }

def save_fasta(sequences: List[Dict], filename: Path, subfamily: str):
    """
    Save sequences in FASTA format.
    
    Args:
        sequences: List of sequence dictionaries
        filename: Output filename
        subfamily: Subfamily identifier
    """
    with open(filename, 'w') as f:
        for seq in sequences:
            # Create descriptive FASTA header
            header = f">{seq['accession']}|{seq['id']}|{subfamily}|{seq['organism']}|{seq['gene']}"
            f.write(f"{header}\n")
            
            # Write sequence in 60-character lines
            sequence = seq['sequence']
            for i in range(0, len(sequence), 60):
                f.write(f"{sequence[i:i+60]}\n")
    
    print(f"  Saved {len(sequences)} sequences to {filename}")

def save_metadata(sequences: List[Dict], filename: Path):
    """
    Save sequence metadata as JSON.
    
    Args:
        sequences: List of sequence dictionaries
        filename: Output filename
    """
    with open(filename, 'w') as f:
        json.dump(sequences, f, indent=2)
    
    print(f"  Saved metadata to {filename}")

def main():
    """Main execution function."""
    print("=" * 80)
    print("Fetching Aromatic Amino Acid Hydroxylase Sequences from UniProt")
    print("=" * 80)
    print()
    
    all_sequences = []
    summary = []
    
    for subfamily_id, config in SUBFAMILIES.items():
        print(f"\n{'─' * 80}")
        print(f"Subfamily: {config['name']} ({subfamily_id})")
        print(f"{'─' * 80}")
        
        # Fetch sequences
        raw_results = fetch_sequences(config['query'], config['max_seqs'])
        
        if not raw_results:
            print(f"  WARNING: No sequences found for {subfamily_id}")
            continue
        
        # Parse results
        sequences = [parse_sequence_record(rec) for rec in raw_results]
        
        # Add subfamily tag
        for seq in sequences:
            seq['subfamily'] = subfamily_id
        
        # Save individual subfamily files
        fasta_file = DATA_DIR / f"{subfamily_id}_sequences.fasta"
        save_fasta(sequences, fasta_file, subfamily_id)
        
        meta_file = DATA_DIR / f"{subfamily_id}_metadata.json"
        save_metadata(sequences, meta_file)
        
        # Collect for combined output
        all_sequences.extend(sequences)
        
        # Summary stats
        avg_length = sum(s['length'] for s in sequences) / len(sequences)
        organisms = set(s['organism'] for s in sequences)
        
        summary.append({
            'subfamily': subfamily_id,
            'name': config['name'],
            'count': len(sequences),
            'avg_length': int(avg_length),
            'organisms': len(organisms)
        })
        
        print(f"  Average length: {avg_length:.0f} aa")
        print(f"  Organisms: {len(organisms)}")
        
        # Be nice to UniProt servers
        time.sleep(1)
    
    # Save combined dataset
    print(f"\n{'=' * 80}")
    print("Saving Combined Dataset")
    print(f"{'=' * 80}")
    
    combined_fasta = DATA_DIR / "all_hydroxylases.fasta"
    save_fasta(all_sequences, combined_fasta, "ALL")
    
    combined_meta = DATA_DIR / "all_hydroxylases_metadata.json"
    save_metadata(all_sequences, combined_meta)
    
    # Print summary
    print(f"\n{'=' * 80}")
    print("SUMMARY")
    print(f"{'=' * 80}")
    print(f"\nTotal sequences retrieved: {len(all_sequences)}")
    print(f"\nBreakdown by subfamily:")
    print(f"{'Subfamily':<10} {'Name':<30} {'Count':<8} {'Avg Length':<12} {'Organisms'}")
    print("─" * 80)
    for s in summary:
        print(f"{s['subfamily']:<10} {s['name']:<30} {s['count']:<8} {s['avg_length']:<12} {s['organisms']}")
    
    print(f"\n{'=' * 80}")
    print("Data retrieval complete!")
    print(f"Files saved to: {DATA_DIR}")
    print(f"{'=' * 80}\n")
    
    # Save summary report
    summary_file = DATA_DIR / "fetch_summary.json"
    with open(summary_file, 'w') as f:
        json.dump({
            'total_sequences': len(all_sequences),
            'subfamilies': summary,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }, f, indent=2)

if __name__ == "__main__":
    main()
