import requests
import json

def get_protein_name_from_pdb(pdb_id):
    """
    Fetches and prints the name of a protein from the RCSB PDB database using its ID.
    """
    # The URL for the RCSB PDB REST API
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    try:
        # Send a GET request to the API
        response = requests.get(url)
        # Raise an exception if the request was unsuccessful
        response.raise_for_status()
        
        # Parse the JSON response
        data = response.json()
        
        # Extract the protein name. The common name is usually stored here.
        # We access the first polymer entity in the structure.
        protein_name = data['rcsb_polymer_entity'][0]['rcsb_polymer_name_combined'][0]['name']
        
        # Extract the full title of the PDB entry for more context
        entry_title = data['struct']['title']
        
        print(f"PDB ID: {pdb_id}")
        print(f"Protein Name: {protein_name}")
        print(f"Full Entry Title: {entry_title}")
        
    except requests.exceptions.RequestException as e:
        print(f"An error occurred with the network request: {e}")
    except (KeyError, IndexError):
        print("Could not parse the protein name from the API response.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# The PDB ID for the protein in the image is '4C2X'
protein_id_from_image = "4C2X"
get_protein_name_from_pdb(protein_id_from_image)