import requests
import json

def get_protein_name(pdb_id):
    """
    Queries the RCSB PDB API to find the name of the protein for a given PDB ID.
    """
    print(f"Querying the Protein Data Bank for ID: {pdb_id}")
    api_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    try:
        response = requests.get(api_url)
        # Raise an exception for bad status codes (4xx or 5xx)
        response.raise_for_status()
        
        data = response.json()
        
        # The official title of the entry usually contains the protein's name
        title = data.get("struct", {}).get("title")
        
        if title:
            print(f"The name of the protein is: {title}")
        else:
            print("Could not find the title for this PDB entry.")

    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err} - Could not find PDB entry {pdb_id}.")
    except Exception as err:
        print(f"An other error occurred: {err}")

# Based on image analysis, the PDB ID for the protein is "1BY8"
protein_id = "1BY8"
get_protein_name(protein_id)