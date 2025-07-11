# First, ensure you have the necessary library installed:
# pip install chembl_webresource_client pandas

from chembl_webresource_client.new_client import new_client
import pandas as pd

def find_interacting_molecules(target_id):
    """
    Finds and prints small molecules that interact with a given ChEMBL target ID.
    """
    print(f"Searching for small molecules that interact with target: {target_id}")
    try:
        # Initialize the client for the 'activity' endpoint
        activity = new_client.activity

        # Query for activities related to the target protein.
        # We filter for activities that have a pChEMBL value, which ensures
        # a quantitative interaction has been measured.
        res = activity.filter(
            target_chembl_id=target_id,
            pchembl_value__isnull=False
        ).only(['molecule_chembl_id'])

        if not res:
            print(f"No interacting small molecules with a pChEMBL value were found for {target_id}.")
            return

        # Convert the list of dictionaries to a pandas DataFrame
        activities_df = pd.DataFrame.from_records(res)

        # Get the unique molecule ChEMBL IDs
        unique_molecule_ids = activities_df['molecule_chembl_id'].unique()

        print(f"\nFound {len(unique_molecule_ids)} unique small molecules interacting with {target_id}.")
        print("Here is a list of their ChEMBL IDs:")
        
        # Print each unique molecule ID
        for molecule_id in unique_molecule_ids:
            print(molecule_id)

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an active internet connection and the 'chembl_webresource_client' library is installed correctly.")

if __name__ == "__main__":
    # The ChEMBL ID for the target protein
    target_chembl_id = 'CHEMBL4803817'
    find_interacting_molecules(target_chembl_id)