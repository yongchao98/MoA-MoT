import requests
import json

def get_protein_info(pdb_id):
    """
    Fetches and prints information about a protein from the RCSB PDB database.

    Args:
        pdb_id (str): The 4-character PDB ID of the protein.
    """
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    print(f"Querying PDB database for ID: {pdb_id}...")
    
    try:
        response = requests.get(url)
        # Raise an exception for bad status codes (4xx or 5xx)
        response.raise_for_status()
        
        data = response.json()
        
        # Extract the title of the structure
        title = data.get('struct', {}).get('title', 'N/A')
        
        # Extract the names of the polymer molecules in the complex
        polymer_entities = data.get('rcsb_polymer_entity', [])
        molecule_names = []
        if polymer_entities:
            for entity in polymer_entities:
                # The names can be in a list, so we join them
                name_list = entity.get('rcsb_polymer_entity_container_identifiers', {}).get('entity_macromolecule_name', [])
                if name_list:
                    molecule_names.append(', '.join(name_list))

        print("\n--- Protein Information ---")
        print(f"PDB ID: {pdb_id}")
        print(f"Structure Title: {title}")
        
        if molecule_names:
            print("\nComponents of the complex:")
            for i, name in enumerate(molecule_names):
                print(f"  Molecule {i+1}: {name}")
        print("-------------------------\n")

    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err} - The PDB ID '{pdb_id}' may not exist.")
    except Exception as err:
        print(f"An error occurred: {err}")

if __name__ == "__main__":
    # The PDB ID identified from the image is 1MOT
    protein_id = "1MOT"
    get_protein_info(protein_id)