# First, you may need to install the ChEMBL web resource client.
# You can do this by running the following command in your terminal:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_interacting_molecules(target_chembl_id):
    """
    Finds and prints small molecules that have recorded interactions
    with a given target ChEMBL ID.
    """
    try:
        # Initialize the client for accessing the ChEMBL database
        activity = new_client.activity

        # Query for activities associated with the specified target ID
        print(f"Querying ChEMBL for activities related to target: {target_chembl_id}...")
        res = activity.filter(target_chembl_id=target_chembl_id)

        if not res:
            print(f"No activity data found for ChEMBL ID: {target_chembl_id}")
            return

        # Use a set to store unique molecule ChEMBL IDs
        molecule_ids = set()
        for item in res:
            # Add the molecule's ChEMBL ID to the set
            if item['molecule_chembl_id']:
                molecule_ids.add(item['molecule_chembl_id'])

        if molecule_ids:
            print(f"\nFound {len(molecule_ids)} unique small molecules that interact with {target_chembl_id}:")
            # Print each molecule ID
            for mol_id in sorted(list(molecule_ids)):
                print(mol_id)
        else:
            print(f"Although activities were found, no specific molecule IDs were associated with them for {target_chembl_id}.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an active internet connection and that the 'chembl_webresource_client' library is installed correctly.")

if __name__ == "__main__":
    target_id = "CHEMBL4803817"
    find_interacting_molecules(target_id)