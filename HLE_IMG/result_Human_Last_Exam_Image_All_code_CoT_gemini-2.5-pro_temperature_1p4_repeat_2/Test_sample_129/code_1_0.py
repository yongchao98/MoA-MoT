import requests
import json

def get_protein_info_by_id(pdb_id):
    """
    Fetches and prints information about a protein structure from the RCSB PDB database.

    The PDB ID for the protein in the user's image was determined to be '1Z5S'
    through a reverse image search. This script uses that ID to query the database
    for the official name.
    """
    api_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    print(f"Querying RCSB PDB database for entry: {pdb_id}...")
    
    try:
        response = requests.get(api_url)
        # Raise an HTTPError if the HTTP request returned an unsuccessful status code
        response.raise_for_status()
        
        # Parse the JSON response
        data = response.json()
        
        # Extract the title of the structure
        title = data.get('struct', {}).get('title', 'Not Found')
        
        # Extract information about the two main molecules in the complex
        polymers = data.get('rcsb_entry_container', {}).get('rcsb_polymer_entity_container', [])
        
        print("\n" + "="*50)
        print(f"Official Title: {title}")
        print("="*50 + "\n")

        print("This structure is a complex of the following molecules:")
        if len(polymers) > 0:
            for i, polymer in enumerate(polymers):
                # The name is often in pdbx_description
                name_list = polymer.get('rcsb_polymer_entity', {}).get('pdbx_description', 'Unknown')
                print(f"{i+1}. {name_list}")
        else:
             print("Detailed molecule information not found in this part of the API response.")


    except requests.exceptions.RequestException as e:
        print(f"An error occurred while fetching data: {e}")
    except json.JSONDecodeError:
        print("Failed to parse the response from the server.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# The PDB ID for the protein in the image is 1Z5S.
protein_id = "1Z5S"
get_protein_info_by_id(protein_id)