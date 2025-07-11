# First, ensure you have the necessary library installed:
# pip install chembl_webresource_client pandas

import pandas as pd
from chembl_webresource_client.new_client import new_client

def find_interacting_molecules(target_id):
    """
    Finds small molecules with measured IC50 values against a given ChEMBL target ID.
    """
    try:
        # Get the activity client from the ChEMBL web service
        activity = new_client.activity
        print(f"Searching for molecules interacting with target: {target_id}")

        # Filter activities for the specific target and standard_type 'IC50'
        # .only() specifies which fields to retrieve to make the query faster
        res = activity.filter(
            target_chembl_id=target_id, 
            standard_type="IC50"
        ).only(['molecule_chembl_id', 'standard_value', 'standard_units'])

        if not res:
            print(f"No molecules with IC50 values found for target {target_id}.")
            return 0
        
        # Convert the list of dictionaries to a pandas DataFrame
        df = pd.DataFrame.from_records(res)

        # Get unique molecule IDs
        unique_molecules = df['molecule_chembl_id'].unique()
        
        count = len(unique_molecules)
        
        print(f"\nFound {count} unique molecules with reported IC50 values against {target_id}:")
        
        # Print each unique molecule ID
        for molecule_id in sorted(list(unique_molecules)):
            print(molecule_id)
        
        return count

    except Exception as e:
        print(f"An error occurred: {e}")
        return 0

if __name__ == "__main__":
    target_chembl_id = "CHEMBL4803817"
    molecule_count = find_interacting_molecules(target_chembl_id)
    # The final answer format is not suitable for a list of molecules.
    # Therefore, I will provide the total count of interacting molecules found.
    # The script itself prints the full list of molecule IDs.
    print(f"\n--- End of List ---")
    print(f"\nFinal Answer: The total number of unique molecules found is {molecule_count}.")
