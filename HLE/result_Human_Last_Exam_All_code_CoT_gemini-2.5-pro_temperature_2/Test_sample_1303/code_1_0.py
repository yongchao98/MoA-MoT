# First, you need to install the required libraries:
# pip install chembl_webresource_client pandas

from chembl_webresource_client.new_client import new_client
import pandas as pd

# The ChEMBL ID for the target (SARS-CoV-2 main protease)
target_chembl_id = 'CHEMBL4803817'

print(f"Searching for small molecule ligands for ChEMBL ID: {target_chembl_id}\n")

try:
    # Initialize the ChEMBL activity client
    activity = new_client.activity

    # Retrieve activities for the target, filtering for:
    # - A valid pChEMBL value (a measure of potency)
    # - Assay type 'B' (for binding assays)
    # We select only the necessary columns to speed up the query.
    res = activity.filter(
        target_chembl_id=target_chembl_id,
        pchembl_value__isnull=False,
        assay_type='B'
    ).only([
        'molecule_chembl_id', 'pchembl_value', 'standard_relation',
        'standard_value', 'standard_units'
    ])

    if not res:
        print(f"No binding activities with pChEMBL values found for target {target_chembl_id}")
    else:
        # Convert the results to a pandas DataFrame for easier manipulation
        df = pd.DataFrame.from_records(res)

        # Ensure 'pchembl_value' is a numeric type for correct sorting
        df['pchembl_value'] = pd.to_numeric(df['pchembl_value'])

        # Sort the DataFrame by pChEMBL value in descending order
        # to find the most potent compounds
        df_sorted = df.sort_values('pchembl_value', ascending=False)

        # Get the top 10 most potent unique ligands
        top_10_ligands = df_sorted.drop_duplicates(subset=['molecule_chembl_id']).head(10)

        # Display the results in a formatted table
        print(f"Top 10 Potent Small Molecules Interacting with {target_chembl_id}:")
        print("-" * 85)
        print(f"{'Molecule ChEMBL ID':<25} | {'pChEMBL Value':<15} | {'Original Measurement (e.g., IC50/Ki)'}")
        print("-" * 85)

        for index, row in top_10_ligands.iterrows():
            # Format the activity string carefully, handling potential None values
            relation = row.get('standard_relation', '') or ''
            value = row.get('standard_value', '') or ''
            units = row.get('standard_units', '') or ''
            activity_info = f"{relation}{value} {units}".strip()

            print(f"{row['molecule_chembl_id']:<25} | {row['pchembl_value']:<15.2f} | {activity_info}")
        
        # Prepare the list of molecule IDs for the final answer
        final_answer_list = top_10_ligands['molecule_chembl_id'].tolist()
        print(f"\n<<<{' , '.join(final_answer_list)}>>>")

except Exception as e:
    print(f"An error occurred: {e}")
    print("Please ensure you have an internet connection and have installed the required libraries:")
    print("pip install chembl_webresource_client pandas")
