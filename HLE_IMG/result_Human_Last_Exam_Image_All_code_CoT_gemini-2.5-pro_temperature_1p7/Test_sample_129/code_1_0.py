import requests
import json

def get_protein_name_from_pdb(pdb_id):
    """
    Retrieves and prints the name of a protein from the RCSB PDB database
    using its PDB ID.
    """
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/A"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)
        
        data = response.json()
        
        # Extract the protein name from the JSON response
        # The 'pdbx_description' field usually contains a descriptive name.
        protein_name = data.get('rcsb_polymer_entity', {}).get('pdbx_description')
        
        if protein_name:
            print(f"The PDB ID for this structure is: {pdb_id}")
            print(f"The name of the protein is: {protein_name}")
        else:
            print(f"Could not find the protein name for PDB ID: {pdb_id}")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while fetching data: {e}")
    except json.JSONDecodeError:
        print("Failed to parse the response as JSON.")
    except KeyError:
        print("Could not find the expected data in the response from the PDB API.")

if __name__ == "__main__":
    # The PDB ID identified for the protein in the image is 1OL5
    pdb_id = "1OL5"
    get_protein_name_from_pdb(pdb_id)