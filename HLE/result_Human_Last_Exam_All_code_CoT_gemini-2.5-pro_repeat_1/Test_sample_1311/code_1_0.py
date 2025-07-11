# First, you may need to install the ChEMBL client library.
# You can do this by running the following command in your terminal:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_binding_affinity():
    """
    Finds the binding affinity (IC50) of SY-1365 to human CDK7 from the ChEMBL database.
    """
    try:
        # Initialize the client
        molecule = new_client.molecule
        target = new_client.target
        activity = new_client.activity

        # --- Step 1: Find the compound (SY-1365) ---
        print("Searching for the compound...")
        # The compound is also known as SY-1365
        compound_name = 'SY-1365'
        molecule_results = molecule.search(compound_name)
        
        if not molecule_results:
            print(f"Could not find the compound '{compound_name}' in the ChEMBL database.")
            return
            
        # Get the ChEMBL ID for the first result
        compound_chembl_id = molecule_results[0]['molecule_chembl_id']
        print(f"Found compound: {compound_name} (ChEMBL ID: {compound_chembl_id})")

        # --- Step 2: Find the target (human CDK7) ---
        print("\nSearching for the target protein...")
        target_name = 'CDK7'
        target_results = target.search(target_name)
        
        # Filter for the human protein
        human_cdk7_target = None
        for t in target_results:
            if t['organism'] == 'Homo sapiens':
                human_cdk7_target = t
                break
        
        if not human_cdk7_target:
            print(f"Could not find the human target '{target_name}' in the ChEMBL database.")
            return

        target_chembl_id = human_cdk7_target['target_chembl_id']
        print(f"Found target: {human_cdk7_target['pref_name']} (ChEMBL ID: {target_chembl_id})")

        # --- Step 3: Retrieve binding affinity data (IC50) ---
        print("\nRetrieving binding affinity data...")
        activity_results = activity.filter(
            molecule_chembl_id=compound_chembl_id,
            target_chembl_id=target_chembl_id,
            standard_type="IC50"  # IC50 is a common measure of potency/affinity
        ).only('standard_type', 'standard_value', 'standard_units', 'assay_description')
        
        if not activity_results:
            print("No IC50 data found for this compound-target pair.")
            return

        # --- Step 4: Display the result ---
        print("\nFound the following binding affinity data:")
        for res in activity_results:
            value = float(res['standard_value'])
            units = res['standard_units']
            print(f"Binding Affinity: IC50 = {value} {units}")
            
            # Determine which range it falls into
            if units == 'nM':
                if 0.1 <= value <= 100:
                    print("This value is in the range of 0.1 - 100 nM.")
                elif value < 0.1:
                    print("This value is in the range of < 0.1 nM.")
            # Add more checks if other units were found
            
    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an internet connection and have installed the 'chembl_webresource_client' library.")

if __name__ == "__main__":
    find_binding_affinity()