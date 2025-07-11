import requests

def get_protein_name_from_pdb(pdb_id):
    """
    Fetches and prints the name of the primary protein entity for a given PDB ID.
    The image provided is a representation of the PDB entry 1L6J.
    """
    api_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    try:
        response = requests.get(api_url)
        # Raise an HTTPError if the HTTP request returned an unsuccessful status code
        response.raise_for_status()
        
        data = response.json()
        
        # The protein's name is located within the entity information.
        # We target the first polymer entity, which corresponds to the main protein chain (MMP-9).
        polymer_entities = data.get('rcsb_entry_container', {}).get('rcsb_polymer_entity', [])
        
        if polymer_entities:
            # The name is stored as the first element in the 'entity_name' list.
            protein_names = polymer_entities[0].get('rcsb_polymer_entity_container_identifiers', {}).get('entity_name', [])
            if protein_names:
                protein_name = protein_names[0]
                print(f"The PDB ID for the structure is: {pdb_id}")
                print(f"The name of the protein is: {protein_name}")
            else:
                print("Protein name could not be found in the entity data.")
        else:
            print("Polymer entity data not found in the API response.")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while communicating with the PDB API: {e}")
    except KeyError:
        print("Failed to parse protein name from the API response due to unexpected data structure.")

# The PDB ID identified from the image is 1L6J
pdb_identifier = "1L6J"
get_protein_name_from_pdb(pdb_identifier)