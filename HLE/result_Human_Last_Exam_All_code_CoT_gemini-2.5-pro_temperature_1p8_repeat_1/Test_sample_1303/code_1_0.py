# First, ensure you have the client library installed:
# pip install chembl_webresource_client pandas

import pandas as pd
from chembl_webresource_client.new_client import new_client

def find_interacting_molecules(target_id):
    """
    Finds and prints small molecules that interact with a given ChEMBL target ID.

    This function queries the ChEMBL database for bioactivities associated with the
    specified target. It filters for common interaction types (IC50, Ki, EC50, Kd)
    and presents the results in a formatted table.

    Args:
        target_id (str): The ChEMBL ID of the target protein.
    """
    print(f"Querying ChEMBL for small molecules interacting with target: {target_id}")
    
    try:
        # Define the resource to query
        activity = new_client.activity

        # Formulate the query to get activities for the target
        # We filter for common standard types to find potent interactions
        # We also filter for a standard_value to ensure data exists
        res = activity.filter(
            target_chembl_id=target_id,
            standard_type__in=['IC50', 'Ki', 'EC50', 'Kd'],
            standard_value__isnull=False
        ).only([
            'molecule_chembl_id', 'standard_type', 'standard_relation', 'standard_value', 'standard_units'
        ])

        if not res:
            print(f"No interactions with specified criteria found for target {target_id}.")
            return

        # Convert the list of dictionaries to a pandas DataFrame for easy handling
        df = pd.DataFrame.from_records(res)

        # Remove duplicate molecules, keeping the first encountered record
        unique_molecules_df = df.drop_duplicates(subset='molecule_chembl_id', keep='first')
        
        # Sort by value for better readability, converting value to numeric
        unique_molecules_df['standard_value'] = pd.to_numeric(unique_molecules_df['standard_value'])
        unique_molecules_df = unique_molecules_df.sort_values(by='standard_value').reset_index(drop=True)


        print("\nFound the following interacting small molecules (showing one entry per molecule):\n")
        # Print the formatted table from the DataFrame
        print(unique_molecules_df[['molecule_chembl_id', 'standard_type', 'standard_relation', 'standard_value', 'standard_units']].to_string())

        # The final answer is the count of unique molecules found
        num_unique_molecules = len(unique_molecules_df)
        print(f"\nTotal number of unique interacting molecules found: {num_unique_molecules}")
        # Final answer format as requested
        print(f"<<<{num_unique_molecules}>>>")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you are connected to the internet and the ChEMBL ID is correct.")


if __name__ == '__main__':
    target_chembl_id = 'CHEMBL4803817'
    find_interacting_molecules(target_chembl_id)