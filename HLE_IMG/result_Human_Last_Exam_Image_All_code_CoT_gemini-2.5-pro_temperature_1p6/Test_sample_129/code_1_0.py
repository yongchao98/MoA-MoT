import requests
import json

def get_protein_name_from_pdb(pdb_id):
    """
    Fetches and prints the name of a protein from the RCSB PDB database
    using its PDB ID.
    """
    # The URL for the PDB REST API endpoint
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"

    try:
        # Make the API request
        response = requests.get(url)
        # Raise an HTTPError if the HTTP request returned an unsuccessful status code
        response.raise_for_status()

        # Parse the JSON response
        data = response.json()

        # Extract the title of the structure
        protein_title = data.get('struct', {}).get('title', 'Title not found')

        print(f"Querying the Protein Data Bank for PDB ID: {pdb_id}")
        print("What is the name of this protein?")
        print(f"The official name for the structure is: {protein_title}")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred during the API request: {e}")
    except json.JSONDecodeError:
        print("Error: Could not decode the JSON response from the server.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # The PDB ID for the protein in the image is 1A5R (Human SUMO-1 protein)
    # The number '1' in '1A5R' is part of the identifier.
    # The numbers in the final output are from the PDB ID itself.
    target_pdb_id = "1A5R"
    get_protein_name_from_pdb(target_pdb_id)