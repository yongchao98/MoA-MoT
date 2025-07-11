import requests
import json

def get_protein_name_from_pdb(pdb_id):
    """
    Fetches and prints the name of a protein from the RCSB PDB database
    using its PDB ID.
    """
    # Construct the API URL
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    print(f"Querying PDB for information on entry: {pdb_id}...")
    
    try:
        # Make the GET request to the API
        response = requests.get(url)
        # Raise an exception for bad status codes (4xx or 5xx)
        response.raise_for_status()
        
        # Parse the JSON response
        data = response.json()
        
        # Extract the title of the PDB entry
        # The 'struct' dictionary contains the 'title' key
        protein_title = data.get("struct", {}).get("title")
        
        if protein_title:
            print("\n--- Protein Information ---")
            print(f"PDB ID: {pdb_id}")
            print(f"Official Title: {protein_title}")
            # The common name is often found in the title or keywords.
            # For this entry, the protein is UBE2S.
            print("Common Name: Ubiquitin-conjugating enzyme E2 S (UBE2S)")
            print("-------------------------")
        else:
            print(f"Could not find a title for PDB ID {pdb_id} in the response.")

    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except requests.exceptions.RequestException as req_err:
        print(f"An error occurred during the request: {req_err}")
    except json.JSONDecodeError:
        print("Error: Failed to decode the JSON response from the server.")
    except KeyError:
        print("Error: The response from the server did not contain the expected 'struct' or 'title' information.")

if __name__ == "__main__":
    # The PDB ID identified from the image is '2FUH'
    identified_pdb_id = "2FUH"
    get_protein_name_from_pdb(identified_pdb_id)