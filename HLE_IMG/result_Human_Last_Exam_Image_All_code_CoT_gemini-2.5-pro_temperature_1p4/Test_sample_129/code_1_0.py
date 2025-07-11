import requests
import json

def get_protein_info(pdb_id):
    """
    Fetches and prints information for a given PDB ID from the RCSB PDB API.
    The protein structure in the image has been identified as PDB entry 1FZY.
    """
    print(f"The protein structure in the image corresponds to the PDB entry: {pdb_id}")
    
    # GraphQL query to get the title and component descriptions
    graphql_query = f"""
    {{
      entry(entry_id: "{pdb_id}") {{
        struct {{
          title
        }}
        polymer_entities {{
          rcsb_polymer_entity {{
            pdbx_description
          }}
        }}
      }}
    }}
    """
    
    api_url = "https://data.rcsb.org/graphql"
    
    try:
        response = requests.post(api_url, json={'query': graphql_query})
        response.raise_for_status()  # Checks for HTTP errors
        data = response.json()
        
        entry_data = data.get('data', {}).get('entry')
        
        if entry_data:
            # Print the official title of the structure
            title = entry_data.get('struct', {}).get('title', 'N/A')
            print(f"\nOfficial Title: {title}")
            
            # Print the names of the protein components
            print("\nComponents:")
            entities = entry_data.get('polymer_entities', [])
            if entities:
                for i, entity in enumerate(entities):
                    description = entity.get('rcsb_polymer_entity', {}).get('pdbx_description', 'N/A')
                    print(f"- Component {i+1}: {description}")
            else:
                print("No polymer components found.")
        else:
            print(f"Could not retrieve detailed information for PDB ID {pdb_id}.")
            
    except requests.exceptions.RequestException as e:
        print(f"An error occurred while contacting the PDB API: {e}")
    except json.JSONDecodeError:
        print("Failed to parse the response from the PDB API.")

# PDB ID identified through visual analysis of the image
identified_pdb_id = "1FZY"
get_protein_info(identified_pdb_id)