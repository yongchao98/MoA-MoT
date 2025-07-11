import requests
import json

def get_protein_name_from_pdb(pdb_id):
    """
    Fetches and prints the name of a protein from the RCSB PDB database using its ID.
    """
    # The URL for the RCSB PDB REST API to get information about a polymer entity
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
    
    try:
        response = requests.get(url)
        # Raise an exception for bad status codes (4xx or 5xx)
        response.raise_for_status()
        
        # Parse the JSON response
        data = response.json()
        
        # Extract the protein name. Based on the API structure, it's under 'rcsb_polymer_entity' -> 'pdbx_description'
        protein_name = data.get('rcsb_polymer_entity', {}).get('pdbx_description')
        
        if protein_name:
            print(f"The protein with PDB ID {pdb_id} is: {protein_name}")
        else:
            print(f"Could not find the name for PDB ID {pdb_id} in the response.")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while fetching data from PDB: {e}")
    except json.JSONDecodeError:
        print("Failed to parse the response from the server.")

# The PDB ID for the protein in the image is 1ZUO
pdb_identifier = "1ZUO"
get_protein_name_from_pdb(pdb_identifier)
