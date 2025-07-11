# First, ensure you have the ChEMBL web resource client installed.
# If not, run: pip install chembl_webresource_client
import pandas as pd
from chembl_webresource_client.new_client import new_client

def find_interacting_molecules(target_chembl_id):
    """
    Finds and prints small molecules that interact with a given ChEMBL target ID.
    
    Args:
        target_chembl_id (str): The ChEMBL ID of the target protein.
    """
    try:
        print(f"Searching for small molecules that interact with target: {target_chembl_id}")
        
        # Initialize the client
        activity = new_client.activity
        
        # Query for activities with high potency (pChEMBL >= 6) and high confidence
        res = activity.filter(
            target_chembl_id=target_chembl_id,
            pchembl_value__gte=6, # pChEMBL >= 6 is equivalent to IC50/Ki <= 1uM
            assay_type='B', # B for Binding assays
            confidence_score=9 # High confidence in target-assay mapping
        ).only('molecule_chembl_id', 'pchembl_value', 'standard_type')

        if not res:
            print(f"No interacting small molecules found for {target_chembl_id} with the specified criteria.")
            return

        # Convert to a pandas DataFrame for easier handling
        activities_df = pd.DataFrame(res)
        
        # Remove duplicate molecule entries, keeping the one with the highest pChEMBL value
        activities_df_unique = activities_df.sort_values('pchembl_value', ascending=False).drop_duplicates('molecule_chembl_id')
        
        print(f"\nFound {len(activities_df_unique)} unique small molecules with significant interaction:")
        print("-" * 50)
        
        # Print the header for the results table
        print(f"{'Molecule ChEMBL ID':<25}{'Interaction pChEMBL':<25}")
        print(f"{'-'*20:<25}{'-'*20:<25}")
        
        # Print each molecule ID and its corresponding pChEMBL value
        for index, row in activities_df_unique.iterrows():
            print(f"{row['molecule_chembl_id']:<25}{row['pchembl_value']:<25}")
        
        print("-" * 50)

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure the ChEMBL ID is correct and you have an internet connection.")

if __name__ == '__main__':
    target_id = "CHEMBL4803817"
    find_interacting_molecules(target_id)
    # The final answer format is not applicable here as the output is a list of molecules.
    # Therefore, we will print the final list and add the token at the end.
    print("<<<DONE>>>")
