import pyensembl
import requests
import json
import os
import pandas as pd
from typing import List, Dict, Any

# ==============================================================================
# Assumptions & Context
# ==============================================================================
# 1. Goal: The assignment asks to identify transcripts and proteins for a specific
#    gene (AKT1, ENSG00000142208) and map the protein coordinates back to the
#    genome.
# 2. Tool Translation: The professor suggested using the R package `ensembldb`.
#    Since you requested a Python solution, this script uses `pyensembl` for local
#    gene/transcript queries and the Ensembl REST API for the more complex
#    protein-to-genome coordinate mapping, which `pyensembl` doesn't handle natively.
# 3. Species & Assembly: Assumed Human (Homo sapiens), Assembly GRCh38 (Release 109),
#    as this is standard unless otherwise specified.
# 4. Input Gene: AKT1 (ENSG00000142208).
# 5. Output: Tab-separated files (TSV) listing transcripts, proteins, and the
#    genomic coordinates for the protein sequences.
# ==============================================================================

# Configuration
GENE_ID = "ENSG00000142208"  # AKT1
GENE_NAME = "AKT1"
ASSEMBLY_RELEASE = 109
SPECIES = "homo_sapiens"

# File output paths
TRANSCRIPTS_OUT = f"{GENE_NAME}_transcripts.tsv"
PROTEINS_OUT = f"{GENE_NAME}_proteins.tsv"
MAPPING_OUT = f"{GENE_NAME}_protein_to_genome_map.tsv"

def setup_pyensembl() -> pyensembl.EnsemblRelease:
    """
    Initializes pyensembl and downloads the necessary reference data if not present.
    """
    print(f"Setting up Ensembl Release {ASSEMBLY_RELEASE}...")
    ensembl = pyensembl.EnsemblRelease(ASSEMBLY_RELEASE)
    
    # Check if we need to download/install data
    try:
        ensembl.gene_by_id(GENE_ID)
    except ValueError:
        print("Reference data not found. Downloading and installing (this may take a while)...")
        ensembl.download()
        ensembl.index()
        print("Installation complete.")
        
    return ensembl

def get_transcripts_and_proteins(ensembl: pyensembl.EnsemblRelease, gene_id: str) -> tuple[List[str], List[str]]:
    """
    Retrieves all transcripts and associated proteins for the given gene.
    """
    print(f"Fetching transcripts and proteins for gene {gene_id}...")
    gene = ensembl.gene_by_id(gene_id)
    transcripts = gene.transcripts
    
    transcript_ids = []
    protein_ids = []
    
    print("\n--- Transcripts & Proteins Found ---")
    for transcript in transcripts:
        t_id = transcript.transcript_id
        p_id = transcript.protein_id
        
        transcript_ids.append(t_id)
        if p_id: # Not all transcripts have a translated protein product
            protein_ids.append(p_id)
            print(f"Transcript: {t_id} -> Protein: {p_id}")
        else:
            print(f"Transcript: {t_id} -> No protein product (non-coding)")
            
    return transcript_ids, protein_ids

def get_protein_to_genome_mapping(protein_id: str) -> List[Dict[str, Any]]:
    """
    Uses the Ensembl REST API to get genomic coordinates for a given protein.
    """
    # Endpoint to map translation (protein) to genome
    server = "https://rest.ensembl.org"
    ext = f"/map/translation/{protein_id}/1..1000000?" # Using a large range to cover the whole protein
    
    url = f"{server}{ext}"
    headers = {"Content-Type": "application/json"}
    
    try:
        response = requests.get(url, headers=headers)
        if not response.ok:
            if response.status_code == 400 and "not a valid" in response.text:
                 # Sometimes obsolete IDs linger, or API limits
                 print(f"  Warning: Protein ID {protein_id} not found via REST API.")
                 return []
            response.raise_for_status()
            
        data = response.json()
        mappings = []
        
        # The API returns mappings of the queried sequence (protein coordinates)
        # to the target sequence (genomic coordinates).
        if 'mappings' in data:
             for m in data['mappings']:
                 # We only care about mapped regions (not gaps/unmapped)
                 if 'mapped' in m:
                    mapped_info = m['mapped']
                    mappings.append({
                        "Protein_ID": protein_id,
                        "Chromosome": mapped_info['seq_region_name'],
                        "Genomic_Start": mapped_info['start'],
                        "Genomic_End": mapped_info['end'],
                        "Strand": mapped_info['strand'],
                        "Protein_Start": m['original']['start'], # AA coordinate start
                        "Protein_End": m['original']['end']      # AA coordinate end
                    })
        return mappings
        
    except requests.exceptions.RequestException as e:
        print(f"  Error fetching mapping for {protein_id}: {e}")
        return []


def main():
    print(f"Starting analysis for {GENE_NAME} ({GENE_ID})")
    
    # 1. Setup PyEnsembl
    ensembl = setup_pyensembl()
    
    # 2. Get Transcripts and Proteins
    transcript_ids, protein_ids = get_transcripts_and_proteins(ensembl, GENE_ID)
    
    # Save Transcripts
    pd.Series(transcript_ids, name="Transcript_ID").to_csv(TRANSCRIPTS_OUT, index=False, header=True)
    print(f"\nSaved {len(transcript_ids)} transcripts to {TRANSCRIPTS_OUT}")
    
    # Save Proteins
    pd.Series(protein_ids, name="Protein_ID").to_csv(PROTEINS_OUT, index=False, header=True)
    print(f"Saved {len(protein_ids)} proteins to {PROTEINS_OUT}")
    
    # 3. Generate Protein to Genome Mapping
    print("\nMapping proteins to genome coordinates using Ensembl REST API...")
    all_mappings = []
    
    for p_id in protein_ids:
        print(f"  Processing {p_id}...")
        mapping_data = get_protein_to_genome_mapping(p_id)
        if mapping_data:
            all_mappings.extend(mapping_data)
            
    if all_mappings:
        df_mappings = pd.DataFrame(all_mappings)
        # Sort for neatness
        df_mappings = df_mappings.sort_values(by=["Protein_ID", "Protein_Start"])
        df_mappings.to_csv(MAPPING_OUT, sep='\t', index=False)
        print(f"\nSaved protein-to-genome mappings to {MAPPING_OUT}")
    else:
        print("\nNo mappings were generated. Please check the API responses.")

    print("\nAnalysis Complete.")
    print("\nNext Step for IGV:")
    print(f"1. Open IGV (installed locally).")
    print(f"2. Select the Human hg38 (GRCh38) genome.")
    print(f"3. Search for {GENE_NAME} in the search bar to view your gene.")

if __name__ == "__main__":
    main()