# First, you may need to install the ChEMBL webresource client.
# You can do this by running the following command in your terminal:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_interacting_molecules(target_chembl_id):
    """
    Finds and prints small molecules that have recorded interactions
    with a given target ChEMBL ID.
    """
    try:
        # Get the activity client
        activity = new_client.activity
        # Search for activities related to the target
        print(f"Querying ChEMBL for molecules interacting with {target_chembl_id}...")
        res = activity.filter(target_chembl_id=target_chembl_id)

        if not res:
            print(f"No interacting molecules found for {target_chembl_id}.")
            return

        # Use a set to store unique molecule IDs
        molecule_ids = set()
        for r in res:
            # We only care about entries that have a molecule ID
            if r['molecule_chembl_id']:
                molecule_ids.add(r['molecule_chembl_id'])

        if molecule_ids:
            print(f"\nFound {len(molecule_ids)} unique small molecules interacting with {target_chembl_id}:")
            # Print each molecule ID
            for mol_id in sorted(list(molecule_ids)):
                print(mol_id)
        else:
            print(f"Found activity entries, but no associated molecule ChEMBL IDs for {target_chembl_id}.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an internet connection and the chembl_webresource_client is installed correctly.")

if __name__ == "__main__":
    target_id = "CHEMBL4803817"
    find_interacting_molecules(target_id)