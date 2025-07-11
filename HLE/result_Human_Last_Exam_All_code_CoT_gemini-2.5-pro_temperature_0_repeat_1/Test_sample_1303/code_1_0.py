# First, ensure you have the client library installed:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client
import pandas as pd

def find_interacting_molecules(target_id):
    """
    Finds and prints small molecules that have recorded interactions
    with a given ChEMBL target ID.
    """
    try:
        # Initialize the client
        activity = new_client.activity

        print(f"Searching for small molecule interactions with target: {target_id}...")

        # Query for activities related to the target, filtering for those with a pChEMBL value
        # A pChEMBL value indicates a quantitative measure of potency
        res = activity.filter(target_chembl_id=target_id, pchembl_value__isnull=False).only(
            ['molecule_chembl_id', 'molecule_pref_name', 'pchembl_value', 'standard_type']
        )

        if not res:
            print(f"No interacting molecules with a pChEMBL value found for {target_id}.")
            return

        # Use a dictionary to store unique molecules by their ID to avoid duplicates
        # and keep their preferred name.
        interacting_molecules = {}
        for r in res:
            # Store the molecule if it's not already there
            if r['molecule_chembl_id'] not in interacting_molecules:
                 interacting_molecules[r['molecule_chembl_id']] = r['molecule_pref_name']

        print(f"\nFound {len(interacting_molecules)} unique small molecules interacting with {target_id}:\n")
        print(f"{'ChEMBL ID':<20} {'Preferred Name'}")
        print("-" * 40)

        # Print each molecule's ID and name
        for chembl_id, name in interacting_molecules.items():
            # If the preferred name is not available, print a placeholder
            display_name = name if name else "N/A"
            print(f"{chembl_id:<20} {display_name}")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an internet connection and the ChEMBL services are available.")
        print("You may also need to install the client library: pip install chembl_webresource_client")


if __name__ == '__main__':
    target_chembl_id = 'CHEMBL4803817'
    find_interacting_molecules(target_chembl_id)