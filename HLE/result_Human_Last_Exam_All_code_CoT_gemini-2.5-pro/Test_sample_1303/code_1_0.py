# First, you may need to install the required libraries:
# pip install chembl_webresource_client pandas

import pandas as pd
from chembl_webresource_client.new_client import new_client

def find_interacting_molecules(target_id, top_n=10):
    """
    Finds small molecules that interact with a given ChEMBL target ID.

    Args:
        target_id (str): The ChEMBL ID of the target protein.
        top_n (int): The number of top molecules to return.

    Returns:
        None. Prints the results.
    """
    print(f"Searching for small molecules interacting with target: {target_id}")
    
    try:
        # Initialize the ChEMBL activity client
        activity = new_client.activity

        # Query for activities related to the target
        # We filter for standard_type 'IC50' and where a pChEMBL value exists
        res = activity.filter(
            target_chembl_id=target_id,
            standard_type="IC50",
            pchembl_value__isnull=False
        ).only([
            'molecule_chembl_id',
            'pchembl_value',
            'standard_value',
            'standard_units'
        ])

        if not res:
            print(f"No molecules with IC50 values found for target {target_id}.")
            return

        # Convert the list of dictionaries to a pandas DataFrame
        df = pd.DataFrame.from_records(res)

        # Data cleaning
        # Ensure pchembl_value is a numeric type for sorting
        df['pchembl_value'] = pd.to_numeric(df['pchembl_value'])
        # Drop rows with non-numeric pchembl_value
        df.dropna(subset=['pchembl_value'], inplace=True)
        # Sort by pChEMBL value in descending order (higher is more potent)
        df_sorted = df.sort_values('pchembl_value', ascending=False)
        # Remove duplicate molecules, keeping the one with the highest potency
        df_unique = df_sorted.drop_duplicates('molecule_chembl_id', keep='first')
        
        # Get the top N results
        top_results = df_unique.head(top_n)

        print(f"\nFound {len(df_unique)} unique molecules.")
        print(f"Displaying the top {len(top_results)} most potent molecules based on pChEMBL value:")
        print("-" * 60)
        print(f"{'Molecule ChEMBL ID':<25} | {'pChEMBL Value':<15} | {'IC50 Value'}")
        print("-" * 60)

        for _, row in top_results.iterrows():
            ic50_val = f"{row['standard_value']} {row['standard_units']}"
            print(f"{row['molecule_chembl_id']:<25} | {float(row['pchembl_value']):<15.2f} | {ic50_val}")
        
        # For the final answer format, we'll return the top molecule
        if not top_results.empty:
            top_molecule_id = top_results.iloc[0]['molecule_chembl_id']
            # This print statement is for the final answer extraction
            print(f"\n<<<Most potent molecule: {top_molecule_id}>>>")


    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an internet connection and the 'chembl_webresource_client' and 'pandas' libraries are installed.")

if __name__ == '__main__':
    target_chembl_id = 'CHEMBL4803817'
    find_interacting_molecules(target_chembl_id)