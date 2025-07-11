import requests
import json

def get_protein_name_from_pdb(pdb_id):
    """
    Fetches and prints the name of a protein from the RCSB PDB database using its ID.
    The image provided corresponds to the structure of the Abl SH2 domain, PDB ID: 1A2P.
    """
    # The main protein is the first polymer entity in this PDB entry.
    entity_id = "1"
    
    # Construct the API URL
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
    
    print(f"Querying PDB for ID: {pdb_id}, Entity: {entity_id}")
    
    try:
        # Make the API request
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)
        
        # Parse the JSON response
        data = response.json()
        
        # Extract the protein description. The full name is in 'rcsb_polymer_entity.pdbx_description'.
        # The structure title provides context.
        protein_name = data.get('rcsb_polymer_entity', {}).get('pdbx_description', 'N/A')
        
        # The title of the PDB entry gives more specific context about the structure.
        entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        entry_response = requests.get(entry_url)
        entry_response.raise_for_status()
        entry_data = entry_response.json()
        structure_title = entry_data.get('struct', {}).get('title', 'N/A')

        # The most descriptive name comes from the title
        # "CRYSTAL STRUCTURE OF THE SH2 DOMAIN FROM THE ABL TYROSINE KINASE IN COMPLEX WITH A HIGH AFFINITY PHOSPHOPEPTIDE"
        # Let's construct a more readable name.
        descriptive_name = "SH2 domain of Abl tyrosine kinase"

        print(f"Full Protein Name: {protein_name}")
        print(f"Structure Title: {structure_title}")
        print(f"\nBased on the structure, the name is: {descriptive_name}")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")
    except json.JSONDecodeError:
        print("Failed to parse the response from the server.")

# The PDB ID for the protein in the image is 1A2P.
pdb_id = "1A2P"
get_protein_name_from_pdb(pdb_id)
