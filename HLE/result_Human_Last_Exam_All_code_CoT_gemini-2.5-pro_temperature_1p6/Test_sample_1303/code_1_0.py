import pandas as pd
from chembl_webresource_client.new_client import new_client

def find_interacting_molecules(target_id: str):
    """
    Finds small molecules that interact with a given ChEMBL target ID.

    Args:
        target_id: The ChEMBL ID of the target protein.
    """
    print(f"Searching for small molecule interactions with target: {target_id}\n")
    try:
        # Initialize the ChEMBL web resource client
        activity = new_client.activity

        # Query for activities related to the target
        # We filter for activities with a standard relation '=' and a pChEMBL value,
        # which indicates a quantified, direct interaction measurement.
        res = activity.filter(target_chembl_id=target_id, pchembl_value__isnull=False, standard_relation='=').only(
            ['molecule_chembl_id', 'type', 'pchembl_value', 'standard_units']
        )

        if not res:
            print(f"No interacting molecules with quantified pChEMBL values found for {target_id}.")
            return

        # Convert the results to a pandas DataFrame
        df = pd.DataFrame(res)

        # Remove duplicates and sort by potency (pChEMBL value)
        df_unique = df.drop_duplicates(subset=['molecule_chembl_id'])
        df_sorted = df_unique.sort_values(by='pchembl_value', ascending=False)
        
        # Limit the output to the top 20 most potent compounds for brevity
        top_results = df_sorted.head(20)

        print(f"Found {len(df_unique)} unique interacting molecules.")
        print("Displaying the top 20 most potent interactions:\n")
        print(f"{'Molecule ChEMBL ID':<20} | {'Activity Type':<15} | {'pChEMBL Value'}")
        print("-" * 60)

        for index, row in top_results.iterrows():
            print(f"{row['molecule_chembl_id']:<20} | {row['type']:<15} | {row['pchembl_value']}")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an internet connection and have installed the necessary library by running:")
        print("pip install chembl_webresource_client pandas")

if __name__ == "__main__":
    target_chembl_id = "CHEMBL4803817"
    find_interacting_molecules(target_chembl_id)