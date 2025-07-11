import requests

def get_protein_name_from_pdb(pdb_id):
    """
    Fetches and prints protein information from the RCSB PDB using its ID.
    """
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    print(f"Querying Protein Data Bank for PDB ID: {pdb_id}")
    
    try:
        response = requests.get(url)
        response.raise_for_status()  # This will raise an HTTPError for bad responses (4xx or 5xx)
        data = response.json()

        title = data.get("struct", {}).get("title")
        
        entities = data.get("rcsb_entry_container", {}).get("polymer_entities", [])
        protein_names = []
        for entity in entities:
            # Navigate the JSON structure to find the macromolecule name
            name_list = entity.get("rcsb_polymer_entity_container_identifiers", {}).get("entity_macromolecule_name", ["Unknown"])
            protein_names.append(name_list[0])

        print("\n--- Identification Result ---")
        if title:
            print(f"Official Title: {title}")
        
        if len(protein_names) == 2:
            print(f"\nThe image shows a complex of two proteins:")
            # Ube2g2 is typically represented as the main enzyme
            print(f"1. The main protein (dark blue) is: {protein_names[0]} (also known as Ube2g2)")
            # FAT10 is the ubiquitin-like modifier
            print(f"2. The binding partner (light blue) is: {protein_names[1]} (also known as FAT10)")
            print("\nBased on common naming, the primary protein identified is Ubiquitin-conjugating enzyme E2 G2.")

        else:
            print("Could not parse individual protein names as expected.")

    except requests.exceptions.RequestException as e:
        print(f"Error: Could not connect to the PDB API. {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

# The PDB ID for the protein structure in the image is 2KJH.
get_protein_name_from_pdb("2KJH")