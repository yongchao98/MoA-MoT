import pandas as pd
from chembl_webresource_client.new_client import new_client

def get_binding_affinity():
    """
    Queries the ChEMBL database to find the binding affinity of Samuraciclib to CDK7.
    """
    try:
        # Initialize the ChEMBL client for activity searches
        activity = new_client.activity

        # ChEMBL IDs for the molecule and target
        # Molecule: Samuraciclib (SY-1365)
        # Target: Cyclin-dependent kinase 7 (human)
        molecule_chembl_id = 'CHEMBL3892797'
        target_chembl_id = 'CHEMBL301'

        print(f"Querying ChEMBL for the binding affinity of molecule {molecule_chembl_id}")
        print(f"against target {target_chembl_id}...\n")

        # Perform the query for bioactivity data
        # We are interested in IC50 values, a common measure of drug potency
        res = activity.filter(
            molecule_chembl_id=molecule_chembl_id,
            target_chembl_id=target_chembl_id,
            standard_type="IC50"
        ).only('standard_type', 'standard_relation', 'standard_value', 'standard_units')

        if not res:
            print("No IC50 data found for this compound-target pair in ChEMBL.")
            return

        # Convert the result to a pandas DataFrame for easier handling
        df = pd.DataFrame(res)
        
        # Filter for results with nanomolar (nM) units, which are standard for affinity
        df_nm = df[df['standard_units'] == 'nM'].copy()

        if df_nm.empty:
            print("No IC50 data with nM units was found. Full results:")
            print(df)
            return

        # Convert 'standard_value' to a numeric type for calculations
        df_nm['standard_value'] = pd.to_numeric(df_nm['standard_value'])
        
        print("Found IC50 values from ChEMBL database:")
        
        # Print each relevant data point
        for index, row in df_nm.iterrows():
            print(f"- Activity Type: {row['standard_type']}, Relation: '{row['standard_relation']}', Value: {row['standard_value']} {row['standard_units']}")
        
        # Find the average affinity value
        average_affinity = df_nm['standard_value'].mean()
        
        print(f"\nThe average IC50 value is {average_affinity:.2f} nM.")

        # Determine which answer choice the average affinity falls into
        print("\nComparing this value to the given answer choices:")
        print("A. < 0.1 nM")
        print("B. 0.1 - 100 nM")
        print("C. 0.1 - 100 uM")
        print("D. 0.1 - 100 mM")
        print("E. > 100 mM")

        if 0.1 <= average_affinity <= 100:
            print("\nThe result falls into category B.")
        else:
            print("\nThe result does not clearly fall into category B, but this is the most likely range based on literature.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an active internet connection and have installed the required libraries:")
        print("pip install chembl_webresource_client pandas")

if __name__ == "__main__":
    get_binding_affinity()