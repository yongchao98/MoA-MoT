import sys
from chembl_webresource_client.new_client import new_client

def find_interacting_molecules(target_id: str):
    """
    Finds and prints small molecules that have recorded bioactivity 
    against a given ChEMBL target ID.
    """
    try:
        # Get the activity resource from the new client
        activity = new_client.activity
        
        # Filter activities by the target ChEMBL ID
        # The client handles pagination automatically
        print(f"Querying ChEMBL for target ID: {target_id}...", file=sys.stderr)
        results = activity.filter(target_chembl_id=target_id)
        
        if not results:
            print(f"No activity data found for ChEMBL ID: {target_id}")
            return
            
        # Use a set to store unique molecule IDs
        molecule_ids = {res['molecule_chembl_id'] for res in results if res['molecule_chembl_id']}
        
        if not molecule_ids:
            print(f"No interacting molecules found for target {target_id}.")
            return

        print(f"\nFound {len(molecule_ids)} unique small molecules interacting with {target_id}:")
        
        # The request asks to print each "number" in the equation.
        # We will interpret this as printing each molecule's ChEMBL ID.
        for mol_id in sorted(list(molecule_ids)):
            print(mol_id)
            
        # The final answer required by the prompt format
        print(f"\n<<<Total unique molecules: {len(molecule_ids)}>>>")


    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        print("Please ensure you have installed the client library with 'pip install chembl_webresource_client'")

if __name__ == "__main__":
    target_chembl_id = "CHEMBL4803817"
    find_interacting_molecules(target_chembl_id)
