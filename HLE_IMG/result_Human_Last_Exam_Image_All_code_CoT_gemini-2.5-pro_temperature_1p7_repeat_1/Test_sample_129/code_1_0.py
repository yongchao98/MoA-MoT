import requests

def get_protein_info(pdb_id):
    """
    Fetches and prints information about a protein structure from the RCSB PDB database.
    """
    print(f"Querying the RCSB PDB database for entry: {pdb_id}\n")
    
    # API endpoint to get the title and entity descriptions
    entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    try:
        # Fetch the main entry data
        response = requests.get(entry_url)
        response.raise_for_status()  # Check for HTTP errors
        entry_data = response.json()

        # Extract and print the title of the structure
        title = entry_data.get('struct', {}).get('title', 'N/A')
        print(f"Official Title: {title}\n")
        
        # Extract and print the names of the individual protein chains (polymers)
        print("The structure contains the following protein components:")
        
        polymer_entities_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}"
        polymer_response = requests.get(f"{polymer_entities_url}/1") # Chain 1
        polymer_response.raise_for_status()
        polymer_data_1 = polymer_response.json()
        description_1 = polymer_data_1.get('rcsb_polymer_entity', {}).get('pdbx_description', 'N/A')
        print(f"- {description_1} (the large, dark blue protein)")

        polymer_response = requests.get(f"{polymer_entities_url}/2") # Chain 2
        polymer_response.raise_for_status()
        polymer_data_2 = polymer_response.json()
        description_2 = polymer_data_2.get('rcsb_polymer_entity', {}).get('pdbx_description', 'N/A')
        print(f"- {description_2} (the small, light blue protein)")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while communicating with the PDB API: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# The PDB ID for the structure in the image is 1W39.
get_protein_info("1W39")