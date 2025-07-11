import requests
import json

def get_protein_info(pdb_id):
    """
    Fetches and prints information about a protein from the RCSB PDB database
    using its PDB ID.
    """
    print(f"--- Querying RCSB PDB for ID: {pdb_id} ---")
    
    # Define the API endpoint for the given PDB ID
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    try:
        # Make the API request
        response = requests.get(url)
        # Raise an error for bad responses (4xx or 5xx)
        response.raise_for_status()
        
        # Parse the JSON data
        data = response.json()
        
        # Extract relevant information from the JSON structure
        struct_title = data.get('struct', {}).get('title', 'Not found')
        protein_name = "Not found"
        
        # The protein name is usually stored in the polymer entity description
        # We'll fetch the first polymer entity associated with the entry
        if data.get('rcsb_entry_container_identifiers', {}).get('polymer_entity_ids'):
            entity_id = data['rcsb_entry_container_identifiers']['polymer_entity_ids'][0]
            entity_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
            entity_response = requests.get(entity_url)
            entity_response.raise_for_status()
            entity_data = entity_response.json()
            protein_name = entity_data.get('rcsb_polymer_entity', {}).get('pdbx_description', 'Not found')
        
        method = data.get('exptl', [{}])[0].get('method', 'Not found')
        
        # Print the retrieved information
        print(f"Protein Name: {protein_name}")
        print(f"Structure Title: {struct_title}")
        print(f"Experimental Method: {method}")
        print("----------------------------------------")
        
        # Conclude with the final name based on the data
        print(f"\nThe protein shown in the image is {protein_name}.")

    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except requests.exceptions.RequestException as err:
        print(f"An error occurred: {err}")
    except json.JSONDecodeError:
        print("Error: Failed to decode the server's response.")
    except KeyError as key_err:
        print(f"Error: Could not find expected key in the data: {key_err}")

# The PDB ID identified from the image is '1NDD'
pdb_id_from_image = "1NDD"
get_protein_info(pdb_id_from_image)