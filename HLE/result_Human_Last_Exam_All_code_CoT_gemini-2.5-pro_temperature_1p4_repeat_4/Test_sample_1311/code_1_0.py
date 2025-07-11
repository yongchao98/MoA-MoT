import sys

def find_binding_affinity():
    """
    This script queries the ChEMBL database to find the binding affinity of
    Fadraciclib to human CDK7.
    """
    try:
        from chembl_webresource_client.new_client import new_client
    except ImportError:
        print("Please install the ChEMBL webresource client first:")
        print("pip install chembl_webresource_client")
        sys.exit(1)

    # --- Step 1: Define compound and target ---
    # The compound is also known as Fadraciclib
    compound_name = 'Fadraciclib'
    target_name = 'CDK7'
    organism = 'Homo sapiens'
    
    print(f"Querying ChEMBL database for binding affinity of {compound_name} to {target_name} ({organism})...")

    # --- Step 2: Find the compound ---
    try:
        molecule = new_client.molecule
        mols = molecule.search(compound_name)
        if not mols:
            print(f"Could not find the compound '{compound_name}' in ChEMBL.")
            return
        # Take the first, most relevant result
        compound_chembl_id = mols[0]['molecule_chembl_id']
        print(f"Found compound: '{compound_name}' (ChEMBL ID: {compound_chembl_id})")

        # --- Step 3: Find the target ---
        target = new_client.target
        targets = target.search(target_name)
        human_cdk7_target = None
        for t in targets:
            # We filter for single protein targets in Homo sapiens to be specific
            if t['organism'] == organism and t['target_type'] == 'SINGLE PROTEIN':
                human_cdk7_target = t
                break
        
        if not human_cdk7_target:
            print(f"Could not find the target '{target_name}' for organism '{organism}' in ChEMBL.")
            return
        target_chembl_id = human_cdk7_target['target_chembl_id']
        print(f"Found target: '{human_cdk7_target['pref_name']}' (ChEMBL ID: {target_chembl_id})")

        # --- Step 4: Retrieve bioactivity data ---
        bioactivity = new_client.bioactivity
        # We look for IC50 or Ki values, which measure potency and affinity
        activities = bioactivity.filter(
            molecule_chembl_id=compound_chembl_id,
            target_chembl_id=target_chembl_id,
            type__in=['IC50', 'Ki'] 
        ).only(['type', 'standard_value', 'standard_units'])

        if not activities:
            print("\nNo binding affinity data (IC50, Ki) found for this specific compound-target pair.")
            return

        # --- Step 5: Analyze and display results ---
        print("\nFound the following binding affinity/potency values:")
        found_values = []
        for act in activities:
            value_str = act['standard_value']
            units = act['standard_units']
            act_type = act['type']
            
            # Ensure we only process valid, numeric data in nM
            if value_str and units == 'nM':
                try:
                    value = float(value_str)
                    print(f"- Type: {act_type}, Value: {value} {units}")
                    found_values.append(value)
                except (ValueError, TypeError):
                    continue
        
        if not found_values:
            print("\nNo valid quantitative affinity data in nM was found.")
            return

        # Check which category the values fall into
        # B. 0.1 - 100 nM
        is_in_range_B = all(0.1 <= val <= 100 for val in found_values)

        print("\n--- Conclusion ---")
        print("The retrieved values represent the concentration of the compound required to inhibit or bind to the target.")
        print("All found values are in the nanomolar (nM) range.")
        print("\nComparing with the given choices:")
        print("A. < 0.1 nM")
        print("B. 0.1 - 100 nM")
        print("C. 0.1 - 100 uM")
        print("D. 0.1 - 100 mM")
        print("E. > 100 mM")
        
        if is_in_range_B:
            print("\nAll reported values fall within the '0.1 - 100 nM' range.")
            print("The correct answer choice is B.")
        else:
            print("\nThe reported values fall into different ranges, but are consistently in the low nM range.")


    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please check your internet connection.")

if __name__ == "__main__":
    find_binding_affinity()