# First, ensure you have the necessary libraries installed.
# You can install them using pip:
# pip install chembl_webresource_client tqdm

from chembl_webresource_client.new_client import new_client
from tqdm import tqdm

def find_interacting_molecules(target_id: str):
    """
    Finds and prints ChEMBL IDs of small molecules that interact with a given target ChEMBL ID.

    Args:
        target_id: The ChEMBL ID of the target protein.
    """
    print(f"Searching for small molecules that interact with target: {target_id}")
    
    try:
        # Initialize the client for the activity resource
        activity = new_client.activity

        # Filter for activities associated with the target that have a pChEMBL value.
        # A pChEMBL value indicates a quantified interaction.
        # This query might return a large number of results, so it can take some time.
        print("Querying ChEMBL database... (this may take a minute)")
        activities = activity.filter(target_chembl_id=target_id, pchembl_value__isnull=False)

        if not activities:
            print(f"No interacting molecules with a pChEMBL value found for {target_id}.")
            return

        # Use a set to store unique molecule IDs
        interacting_molecule_ids = set()

        # Use tqdm for a progress bar as we iterate through potentially many results
        for act in tqdm(activities, desc="Processing activities"):
            if 'molecule_chembl_id' in act:
                interacting_molecule_ids.add(act['molecule_chembl_id'])

        if interacting_molecule_ids:
            print(f"\nFound {len(interacting_molecule_ids)} unique interacting small molecules:")
            # Sort the IDs for consistent output and print each one
            for molecule_id in sorted(list(interacting_molecule_ids)):
                print(molecule_id)
        else:
            print("Found activities, but could not extract any molecule ChEMBL IDs.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an active internet connection and that the ChEMBL ID is correct.")

if __name__ == "__main__":
    # The ChEMBL ID for the target protein (SARS-CoV-2 Main protease)
    target_chembl_id = "CHEMBL4803817"
    find_interacting_molecules(target_chembl_id)