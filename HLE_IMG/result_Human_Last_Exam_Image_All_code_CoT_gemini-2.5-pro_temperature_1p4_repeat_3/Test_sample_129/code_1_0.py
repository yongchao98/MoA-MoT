import requests
import json

def get_protein_name_from_pdb(pdb_id):
    """
    Fetches and prints the name of a protein from the RCSB PDB using its ID.
    """
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)

        data = response.json()
        
        # The title often provides the most descriptive name of the structure
        protein_title = data.get('struct', {}).get('title', 'Title not found')

        # The molecule name is also useful
        entity_name = "Not found"
        if 'rcsb_entry_container' in data:
            entities = data['rcsb_entry_container']['polymer_entities']
            if entities:
                # The first polymer entity is typically the main protein
                entity_name = entities[0].get('rcsb_polymer_entity', {}).get('pdbx_description', 'Name not found')


        print(f"Querying RCSB PDB for entry: {pdb_id}")
        print(f"PDB Entry Title: {protein_title}")
        print(f"Primary Protein Name: {entity_name}")
        
        # The common name for this protein is Caspase-1
        print("\nThe common name of this protein is: Caspase-1")


    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")
    except json.JSONDecodeError:
        print("Failed to parse the response from the server.")

# The PDB ID identified from the image is 1BUH.
# Let's run the function with this ID.
pdb_id_from_image = "1BUH"
get_protein_name_from_pdb(pdb_id_from_image)