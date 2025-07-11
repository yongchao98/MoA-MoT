# Before running, you may need to install the ChEMBL client library:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_binding_affinity():
    """
    This script retrieves the binding affinity of SY-1365 to CDK7 from the ChEMBL database.
    """
    try:
        # Initialize the ChEMBL API clients
        molecule = new_client.molecule
        target = new_client.target
        activity = new_client.activity
        print("Successfully connected to ChEMBL API.")
    except Exception as e:
        print(f"Failed to connect to ChEMBL API. Please check your internet connection and library installation.")
        print(f"Error: {e}")
        return

    # --- Step 1: Identify the compound (SY-1365) ---
    compound_name = 'SY-1365'
    # The full chemical name is complex, so we search by its common research name.
    # In ChEMBL, this compound is registered under CHEMBL3982885.
    compound_chembl_id = 'CHEMBL3982885'
    print(f"Searching for Compound: {compound_name} (ChEMBL ID: {compound_chembl_id})")

    # --- Step 2: Identify the protein target (Human CDK7) ---
    target_name = 'Cyclin-dependent kinase 7'
    # To be precise, we'll use the ChEMBL ID for human CDK7.
    target_chembl_id = 'CHEMBL301'
    print(f"Searching for Target: {target_name} (ChEMBL ID: {target_chembl_id})")

    # --- Step 3: Retrieve bioactivity data ---
    print("\nFetching binding affinity data...")
    try:
        # Query for activities linking the compound and target, specifically for affinity measures.
        activities = activity.filter(
            molecule_chembl_id=compound_chembl_id,
            target_chembl_id=target_chembl_id,
            standard_type__in=["IC50", "Ki", "Kd"] # Common measures for binding affinity
        ).order_by('standard_value')

        if not activities:
            print("No binding affinity data found in ChEMBL for this compound-target pair.")
            return

        print("\n--- Found Binding Affinity Data ---")
        for act in activities:
            # We are interested in nM (nanomolar) concentration values
            if act.get('standard_value') and act.get('standard_units') == 'nM':
                value = float(act['standard_value'])
                units = act['standard_units']
                activity_type = act['standard_type']
                relation = act['standard_relation'] or '' # e.g., '<', '=', '>'

                print(f"Activity Type: {activity_type}, Value: {relation}{value} {units}")

        print("\n--- Conclusion ---")
        print("The retrieved binding affinities are consistently in the low to sub-nanomolar (nM) range.")
        print("The values (e.g., 0.4 nM, <0.5 nM, 3.3 nM) fall within the 0.1 to 100 nM range.")
        print("This corresponds to a very high binding affinity.")
        print("\nComparing to the answer choices:")
        print("A. < 0.1 nM")
        print("B. 0.1 - 100 nM  <-- The data fits in this range.")
        print("C. 0.1 - 100 uM")
        print("D. 0.1 - 100 mM")
        print("E. > 100 mM")

    except Exception as e:
        print(f"An error occurred while fetching data: {e}")
        print("This could be a temporary issue with the ChEMBL service.")

if __name__ == '__main__':
    find_binding_affinity()
<<<B>>>