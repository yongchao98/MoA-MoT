# Before running, make sure you have the necessary libraries installed:
# pip install chembl-webresource-client pandas

import pandas as pd
from chembl_webresource_client.new_client import new_client

def find_interacting_molecules(target_chembl_id: str):
    """
    Finds and displays small molecules that interact with a given ChEMBL target ID.
    """
    print(f"Searching for small molecules interacting with target: {target_chembl_id}\n")

    try:
        # Initialize the client for the activity resource
        activity = new_client.activity

        # Filter for activities related to the target ID that have a pChEMBL value
        # pChEMBL is a standardized measure of potency (-log10(IC50, Ki, etc.))
        res = activity.filter(
            target_chembl_id=target_chembl_id
        ).filter(
            pchembl_value__isnull=False
        )

        if not res:
            print(f"No activities with a pChEMBL value found for {target_chembl_id}.")
            return None

        # Convert the results to a pandas DataFrame
        df = pd.DataFrame.from_records(res)

        # Select and rename important columns for clarity
        df_selected = df[[
            'molecule_chembl_id',
            'canonical_smiles',
            'standard_type',
            'pchembl_value'
        ]]

        # Convert pchembl_value to numeric and sort by it
        df_selected['pchembl_value'] = pd.to_numeric(df_selected['pchembl_value'])
        df_sorted = df_selected.sort_values(by='pchembl_value', ascending=False).reset_index(drop=True)

        print(f"Found {len(df_sorted)} interacting molecules.")
        print("Displaying the Top 20 most potent molecules:")
        print(df_sorted.head(20).to_string())

        # The most potent molecule is the first one in the sorted list
        most_potent_molecule_id = df_sorted.iloc[0]['molecule_chembl_id']
        return most_potent_molecule_id

    except Exception as e:
        print(f"An error occurred: {e}")
        return None

if __name__ == '__main__':
    target_id = "CHEMBL4803817"
    most_potent_id = find_interacting_molecules(target_id)
    if most_potent_id:
        print(f"\nThe most potent interacting molecule is: {most_potent_id}")
        # Final answer format
        print(f"\n<<<{most_potent_id}>>>")
