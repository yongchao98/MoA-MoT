import requests
import json

def get_protein_info(pdb_id):
    """
    Fetches and prints information about a protein structure from the RCSB PDB database.
    
    Args:
        pdb_id (str): The 4-character PDB ID.
    """
    print(f"--- Querying RCSB PDB for entry {pdb_id} ---")
    
    # URL for the PDB API to get descriptive information about an entry
    summary_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    try:
        # Fetch the main entry data
        response = requests.get(summary_url)
        response.raise_for_status()  # Raise an error for bad responses (4xx or 5xx)
        data = response.json()
        
        # Print the title of the experiment
        title = data.get('struct', {}).get('title', 'Title not found.')
        print(f"\nPDB Entry Title: {title}\n")
        
        # Get the list of polymer entity IDs
        polymer_entities = data.get('rcsb_entry_container_identifiers', {}).get('polymer_entity_ids', [])
        
        if not polymer_entities:
            print("No polymer entities found in this entry.")
            return
            
        print("Protein components identified in this structure:")
        # Loop through each polymer entity to get its name
        for entity_id in polymer_entities:
            # Construct URL to get specific polymer entity information
            entity_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
            entity_response = requests.get(entity_url)
            entity_response.raise_for_status()
            entity_data = entity_response.json()
            
            # Extract the descriptive name of the molecule
            molecule_name = entity_data.get('rcsb_polymer_entity', {}).get('pdbx_description', 'Name not available')
            print(f"- Molecule {entity_id}: {molecule_name}")

    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err} - Could not retrieve data for PDB ID '{pdb_id}'.")
    except Exception as err:
        print(f"An error occurred: {err}")

# The PDB ID for the protein in the image is 1NW9
protein_pdb_id = "1NW9"
get_protein_info(protein_pdb_id)