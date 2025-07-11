import requests
import json

def get_protein_info(pdb_id):
    """
    Fetches and prints protein information from the RCSB PDB API.
    
    Args:
        pdb_id (str): The 4-character PDB ID.
    """
    # For this entry, there is one protein chain (polymer entity), with entity_id '1'
    entity_id = "1"
    api_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
    
    print(f"Querying RCSB PDB API for entry: {pdb_id}...")
    
    try:
        response = requests.get(api_url)
        # Raise an HTTPError for bad responses (4xx or 5xx)
        response.raise_for_status()
        
        data = response.json()
        
        # Extract the common name/description of the protein
        protein_name = data.get("rcsb_polymer_entity", {}).get("pdbx_description", "Name not found")
        
        # Extract the source organism
        organism_list = data.get("rcsb_entity_source_organism", [])
        if organism_list:
            organism = organism_list[0].get("ncbi_scientific_name", "Organism not found")
        else:
            organism = "Organism not found"
            
        print("\n--- Protein Information ---")
        print(f"Name: {protein_name}")
        print(f"Source Organism: {organism}")
        print(f"PDB ID: {pdb_id}")
        print("-------------------------")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while connecting to the PDB API: {e}")
    except json.JSONDecodeError:
        print("Error: Failed to decode the response from the API.")
    except (KeyError, IndexError):
        print("Error: Could not find the expected information in the API response.")

if __name__ == "__main__":
    # The PDB ID for the protein in the image is 1J7D
    protein_pdb_id = "1J7D"
    get_protein_info(protein_pdb_id)