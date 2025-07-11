import requests

def get_protein_name_from_pdb(pdb_id):
    """
    Fetches and prints the name of a protein from the RCSB PDB database
    using its PDB ID.

    Args:
        pdb_id (str): The 4-character PDB ID.
    """
    # We query the 'polymer_entity' endpoint for the first entity (ID '1')
    # to get the specific name of the protein molecule.
    api_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
    
    print(f"Querying RCSB PDB API for PDB ID: {pdb_id}...")

    try:
        response = requests.get(api_url)
        # Raise an error for bad responses (e.g., 404 Not Found, 500 Server Error)
        response.raise_for_status()
        
        data = response.json()
        
        # The protein's descriptive name is located in the 'pdbx_description' field
        # within the 'rcsb_polymer_entity' object.
        protein_name = data.get('rcsb_polymer_entity', {}).get('pdbx_description')
        
        if protein_name:
            print(f"Successfully retrieved protein name.")
            print(f"The name of the protein is: {protein_name}")
        else:
            print("Could not find the protein name in the API response.")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while communicating with the API: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# The PDB ID identified from the image is '1A2P'.
protein_id = "1A2P"
get_protein_name_from_pdb(protein_id)