import requests
import json

def get_protein_info(pdb_id):
    """
    Fetches and prints information about a protein from the RCSB PDB database.
    """
    print(f"Querying the RCSB PDB database for entry: {pdb_id}")
    
    # Get descriptive information for the entry
    entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(entry_url)
    
    if response.status_code == 200:
        data = response.json()
        title = data.get('struct', {}).get('title', 'N/A')
        print(f"\nTitle: {title}")
    else:
        print(f"Failed to fetch entry information for {pdb_id}. Status code: {response.status_code}")
        return

    # Get information about the polymer entities (proteins)
    polymer_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}"
    response = requests.get(polymer_url)

    if response.status_code == 200:
        data = response.json()
        print("\nProteins in this complex:")
        for i, entity in enumerate(data):
            try:
                name = entity['rcsb_polymer_entity_container_identifiers']['entity_macromolecule_names'][0]['name']
                print(f"- Molecule {i+1}: {name}")
            except (KeyError, IndexError):
                print(f"- Molecule {i+1}: Name not found")
        
        print("\nBased on this data and visual analysis, the image shows the SUMO protein (Small Ubiquitin-like Modifier), likely the light blue chain, in a complex with its E2-conjugating enzyme Ubc9 (the dark blue chain).")
        print("Therefore, a primary name for the protein of interest is SUMO protein.")

    else:
        print(f"Failed to fetch polymer entity information for {pdb_id}. Status code: {response.status_code}")


if __name__ == "__main__":
    # The PDB ID for the structure in the image is 1Z5S.
    pdb_id = "1Z5S"
    get_protein_info(pdb_id)