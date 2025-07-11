# First, ensure you have the client installed:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_binding_affinity():
    """
    Finds the binding affinity of a compound to a target protein
    by querying the ChEMBL database.
    """
    try:
        # Initialize the ChEMBL client for activities
        activity = new_client.activity

        # --- Step 1: Search for the compound by its synonym 'SY-1365' ---
        molecule_search = new_client.molecule.search('SY-1365')
        if not molecule_search:
            print("Could not find the compound SY-1365 in the ChEMBL database.")
            return
        
        # Take the first and most likely result
        compound = molecule_search[0]
        compound_chembl_id = compound['molecule_chembl_id']
        print(f"Found Compound: SY-1365 (ChEMBL ID: {compound_chembl_id})")

        # --- Step 2: Search for the target 'CDK7' for Homo sapiens ---
        target_search = new_client.target.filter(pref_name__iexact='Cyclin-dependent kinase 7', organism='Homo sapiens')
        if not target_search:
            print("Could not find the target CDK7 (Homo sapiens) in the ChEMBL database.")
            return

        target = target_search[0]
        target_chembl_id = target['target_chembl_id']
        print(f"Found Target: {target['pref_name']} (ChEMBL ID: {target_chembl_id})")

        # --- Step 3: Query for binding activity data (IC50) ---
        print("\nSearching for binding affinity data...")
        
        # Filter for IC50 values measured in binding assays ('B')
        res = activity.filter(
            target_chembl_id=target_chembl_id,
            molecule_chembl_id=compound_chembl_id,
            standard_type="IC50",
            assay_type='B'
        ).only(['standard_type', 'standard_value', 'standard_units', 'standard_relation'])

        if not res:
            print(f"No IC50 binding data found for {compound_chembl_id} on {target_chembl_id}.")
            return
            
        # --- Step 4: Print the first relevant result ---
        for entry in res:
            if entry['standard_value'] and entry['standard_units'] == 'nM':
                affinity_value = float(entry['standard_value'])
                relation = entry['standard_relation']
                print("\n--- Binding Affinity Result ---")
                print(f"The binding affinity (IC50) is: {relation} {affinity_value} {entry['standard_units']}")
                # This value of ~6 nM falls within the 0.1 - 100 nM range.
                return

        print("No IC50 data with 'nM' units found.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an active internet connection and the 'chembl_webresource_client' library is installed ('pip install chembl_webresource_client').")

if __name__ == "__main__":
    find_binding_affinity()