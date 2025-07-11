import requests
import json

def get_protein_name_from_pdb(pdb_id):
    """
    Fetches the protein name from the RCSB PDB using its ID.
    """
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    try:
        response = requests.get(url)
        # Raise an exception for bad status codes (4xx or 5xx)
        response.raise_for_status()
        
        data = response.json()
        
        # The title of the PDB entry describes the structure
        protein_name = data.get("struct", {}).get("title")
        
        if protein_name:
            print(f"The PDB entry {pdb_id} corresponds to the following structure:")
            print(protein_name)
        else:
            print(f"Could not find the name for PDB ID: {pdb_id}")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while fetching data from the PDB: {e}")
    except json.JSONDecodeError:
        print("Failed to parse the response from the PDB API.")

if __name__ == "__main__":
    # The PDB ID for the protein in the image is 1k9c.
    pdb_id = "1k9c"
    get_protein_name_from_pdb(pdb_id)