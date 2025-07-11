# First, ensure you have the necessary libraries installed.
# If not, you can install them by running:
# pip install chembl_webresource_client pandas

from chembl_webresource_client.new_client import new_client
import pandas as pd

def find_interacting_molecules(target_id):
    """
    Finds and prints small molecules that interact with a given ChEMBL target ID.
    """
    print(f"Querying ChEMBL for target: {target_id}")

    try:
        # 1. Initialize the ChEMBL activity client
        activity = new_client.activity

        # 2. Query for activities with a pChEMBL value for the specified target
        res = activity.filter(
            target_chembl_id=target_id,
            pchembl_value__isnull=False  # Ensures we get activities with a potency value
        )

        if not res:
            print(f"No activities with pChEMBL values found for target {target_id}.")
            return 0

        # 3. Convert the results to a pandas DataFrame for easier processing
        df = pd.DataFrame.from_records(res)

        # 4. Extract the unique molecule ChEMBL IDs
        unique_molecules = df['molecule_chembl_id'].unique()
        count = len(unique_molecules)

        if count > 0:
            print(f"\nFound {count} unique small molecules interacting with {target_id}:")
            # 5. Print each unique molecule ID
            for molecule_id in unique_molecules:
                print(molecule_id)
        else:
            print("No unique interacting molecules found.")

        return count

    except Exception as e:
        print(f"An error occurred: {e}")
        return 0

if __name__ == '__main__':
    target_chembl_id = 'CHEMBL4803817'
    total_count = find_interacting_molecules(target_chembl_id)

    # Final answer format
    print(f"\nTotal number of unique interacting molecules: {total_count}")
    print(f"<<<{total_count}>>>")
