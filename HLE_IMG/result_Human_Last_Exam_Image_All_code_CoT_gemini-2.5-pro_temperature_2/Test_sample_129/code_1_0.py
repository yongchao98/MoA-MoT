import requests
import json

def get_protein_info(pdb_id):
    """
    Fetches and prints information for a given PDB ID from the RCSB PDB database.
    """
    print(f"Querying the RCSB PDB database for PDB ID: {pdb_id}")
    
    # RCSB PDB API URL for entry information
    entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    # RCSB PDB API URL for polymer entity information
    polymer_entity_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1" # For chain A (Cbl-b)
    ubiquitin_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/2" # For chain B (Ubiquitin)
    
    try:
        # Fetching general structure information
        response_entry = requests.get(entry_url)
        response_entry.raise_for_status()  # Raise an exception for bad status codes
        data_entry = response_entry.json()
        
        # Fetching info for the first polymer (Cbl-b TKB domain)
        response_poly1 = requests.get(polymer_entity_url)
        response_poly1.raise_for_status()
        data_poly1 = response_poly1.json()

        # Fetching info for the second polymer (Ubiquitin)
        response_poly2 = requests.get(ubiquitin_url)
        response_poly2.raise_for_status()
        data_poly2 = response_poly2.json()

        print("\n--- Protein Information ---")
        
        # Print the official title of the structure
        title = data_entry.get('struct', {}).get('title', 'N/A')
        print(f"Structure Title: {title}")
        
        # Print information about the main protein components
        protein1_name = data_poly1.get('rcsb_polymer_entity', {}).get('pdbx_description', 'N/A')
        protein2_name = data_poly2.get('rcsb_polymer_entity', {}).get('pdbx_description', 'N/A')

        print(f"\nThe structure is a complex containing two main components:")
        print(f"1. Component 1 (dark blue in image): {protein1_name}")
        print(f"2. Component 2 (light blue in image): {protein2_name}")
        
        # A more common name for E3 ubiquitin-protein ligase CBL-B is Cbl-b
        print("\nCommon Name: The protein shown is the TKB domain of Cbl-b covalently attached to Ubiquitin.")
        
    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except Exception as err:
        print(f"An error occurred: {err}")

if __name__ == '__main__':
    # The PDB ID for the structure in the image is 1ZBK
    pdb_id = "1ZBK"
    get_protein_info(pdb_id)