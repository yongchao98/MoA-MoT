import requests
import json

def get_protein_name_by_id(pdb_id):
    """
    Fetches the full name of a protein from the RCSB PDB database using its ID
    and prints it.
    """
    # Construct the API URL for the specified PDB entry
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    print(f"Querying RCSB PDB database for entry: {pdb_id}...")
    
    try:
        # Send a GET request to the PDB API
        response = requests.get(url)
        # Raise an exception for bad status codes (4xx or 5xx)
        response.raise_for_status()
        
        # Parse the JSON response
        data = response.json()
        
        # Extract the title from the structure information
        # The full name is typically stored in the 'struct.title' field
        protein_name = data.get('struct', {}).get('title', 'Name not found.')
        
        print("\n--- Protein Information ---")
        print(f"PDB ID: {pdb_id}")
        print(f"Name: {protein_name}")
        print("---------------------------")
        
    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except requests.exceptions.RequestException as err:
        print(f"An error occurred: {err}")
    except json.JSONDecodeError:
        print("Error: Failed to parse the response from the server.")

# The PDB ID for the protein in the image is 1CMW
protein_pdb_id = "1CMW"
get_protein_name_by_id(protein_pdb_id)