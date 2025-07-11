# First, ensure you have the ChEMBL web resource client installed:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_interacting_molecules(target_id):
    """
    Finds and prints small molecules with reported IC50 activity
    against a specific ChEMBL target ID.
    """
    print(f"Querying ChEMBL for small molecules interacting with target: {target_id}\n")
    
    try:
        # Get the activity client
        activity = new_client.activity

        # Filter activities for the target
        # We look for small molecules with a reported IC50 value
        activities = activity.filter(
            target_chembl_id=target_id,
            standard_type="IC50",
            molecule_type="Small molecule"
        )

        if not activities:
            print(f"No small molecules with reported IC50 values found for target {target_id}.")
            return

        # Use a set to store unique molecule IDs to avoid duplicates
        molecule_ids = set()
        for act in activities:
            molecule_ids.add(act['molecule_chembl_id'])

        # Sort the IDs for consistent output
        sorted_molecule_ids = sorted(list(molecule_ids))

        print(f"Found {len(sorted_molecule_ids)} unique small molecules interacting with {target_id}:")
        for mol_id in sorted_molecule_ids:
            print(mol_id)

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure the chembl_webresource_client is installed ('pip install chembl_webresource_client') and you have an internet connection.")

if __name__ == "__main__":
    target_chembl_id = "CHEMBL4803817"
    find_interacting_molecules(target_chembl_id)