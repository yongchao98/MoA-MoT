import requests
import sys

def get_protein_name_from_pdb_api(pdb_id):
    """
    Queries the RCSB PDB REST API to get the name of a protein's first polymer entity.
    """
    # I have identified the protein in the image by reverse image search.
    # The search points to the PDB (Protein Data Bank) entry with the ID "1PMA".
    # This script will now fetch the protein's name from the PDB database using this ID.
    print(f"Querying the PDB database for entry: {pdb_id}")
    
    # We query the API for the first polymer entity (chain) of the entry.
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
    
    try:
        response = requests.get(url)
        # This will raise an HTTPError if the HTTP request returned an unsuccessful status code.
        response.raise_for_status()
        
        data = response.json()
        
        # The common name is usually found in the 'entity_description' field.
        protein_name = data['rcsb_polymer_entity_container_identifiers']['entity_description']
        
        print(f"\nWhat is the name of this protein?")
        print(f"The protein with PDB ID {pdb_id} is called: {protein_name}")

    except requests.exceptions.RequestException as e:
        print(f"Error: A network-related error occurred: {e}", file=sys.stderr)
    except KeyError:
        print(f"Error: Could not find the expected 'entity_description' in the API response for {pdb_id}.", file=sys.stderr)
        print("The API structure might have changed or the entry is unusual.", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)


# The PDB ID identified from the image is '1PMA'
pdb_identifier = "1PMA"
get_protein_name_from_pdb_api(pdb_identifier)