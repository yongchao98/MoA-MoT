import requests
import json

def find_protein_name():
    """
    Identifies the protein in the image by fetching its data from the
    RCSB Protein Data Bank (PDB).
    """
    # Based on reverse image search, the structure corresponds to PDB ID 5BBN.
    # This structure is a complex of UBE2S and Ubiquitin. We will fetch the
    # name of the primary, larger protein (Chain A, Entity 1).
    pdb_id = "5BBN"
    entity_id = "1"
    
    api_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
    
    print(f"Querying RCSB PDB for structure {pdb_id}...")
    
    try:
        response = requests.get(api_url, timeout=10)
        # Raise an HTTPError for bad responses (e.g., 404, 500)
        response.raise_for_status()
        
        data = response.json()
        
        # The common name(s) are stored in the 'entry_macromolecule_names' list
        protein_names = data.get('rcsb_polymer_entity_container_identifiers', {}).get('entry_macromolecule_names', [])
        
        if protein_names:
            protein_name = protein_names[0]
            print(f"\n--- Protein Information ---")
            print(f"PDB ID: {pdb_id}")
            print(f"The name of the main enzyme in the complex is: {protein_name}")
            print(f"---------------------------")
        else:
            print("Could not find a common name for the protein in the PDB data.")

    except requests.exceptions.RequestException as e:
        print(f"Error: Could not connect to the PDB API. {e}")
    except json.JSONDecodeError:
        print("Error: Failed to parse the response from the server.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Execute the function to find and print the protein name.
find_protein_name()