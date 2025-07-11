import json
import urllib.request

def get_protein_name_from_pdb(pdb_id):
    """
    Fetches protein information from the RCSB PDB database and returns its name.

    Args:
        pdb_id (str): The 4-character PDB accession code.

    Returns:
        str: The name of the protein, or an error message if not found.
    """
    # The REST API URL to get core entry data for a given PDB ID
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    try:
        # Send a request to the API
        with urllib.request.urlopen(url) as response:
            if response.status == 200:
                # Read and parse the JSON response
                data = json.loads(response.read().decode('utf-8'))
                
                # The protein name is located in the polymer_entities section
                # We are targeting the first polymer entity in the list
                protein_name = data['polymer_entities'][0]['rcsb_polymer_entity']['pdbx_description']
                
                print(f"The protein corresponding to PDB ID {pdb_id} is: {protein_name}")
            else:
                print(f"Error: Received status code {response.status} for PDB ID {pdb_id}")

    except urllib.error.HTTPError as e:
        print(f"Error fetching data for PDB ID {pdb_id}: {e}")
    except (KeyError, IndexError):
        print(f"Could not parse the protein name from the response for PDB ID {pdb_id}.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# The PDB ID identified from the image is 1P4O
protein_id = "1P4O"
get_protein_name_from_pdb(protein_id)