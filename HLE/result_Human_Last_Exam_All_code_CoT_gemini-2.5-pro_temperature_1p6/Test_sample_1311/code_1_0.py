# First, ensure you have the necessary library installed:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_binding_affinity():
    """
    Finds the binding affinity (IC50) of a compound (SY-1365) to a target (CDK7)
    by querying the ChEMBL database.
    """
    try:
        # Initialize clients for different ChEMBL resources
        molecule = new_client.molecule
        target = new_client.target
        activity = new_client.activity

        print("Step 1: Searching for the compound SY-1365...")
        # Find the molecule using its preferred name (synonym)
        mol_results = molecule.filter(pref_name__iexact='SY-1365').only('molecule_chembl_id', 'pref_name')
        if not mol_results:
            print("Compound SY-1365 not found.")
            return

        compound = mol_results[0]
        mol_chembl_id = compound['molecule_chembl_id']
        print(f"Found Compound: {compound['pref_name']} (ID: {mol_chembl_id})")

        print("\nStep 2: Searching for the target Human CDK7...")
        # Find the human target protein
        target_results = target.filter(pref_name__iexact='Cyclin-dependent kinase 7', organism='Homo sapiens').only('target_chembl_id', 'pref_name')
        if not target_results:
            print("Target CDK7 (Human) not found.")
            return

        protein_target = target_results[0]
        target_chembl_id = protein_target['target_chembl_id']
        print(f"Found Target: {protein_target['pref_name']} (ID: {target_chembl_id})")

        print("\nStep 3: Querying for binding affinity (IC50)...")
        # Find activity data linking the compound and target, filtering for IC50 values
        activity_results = activity.filter(
            molecule_chembl_id=mol_chembl_id,
            target_chembl_id=target_chembl_id,
            standard_type="IC50"
        ).only('standard_type', 'standard_value', 'standard_units', 'relation')

        if not activity_results:
            print("No IC50 binding affinity data found for this pair.")
            return

        # Print the first relevant result found
        result = activity_results[0]
        print("\n--- Final Result ---")
        print(f"Binding Affinity Found for {compound['pref_name']} to {protein_target['pref_name']}:")
        print(f"Measurement Type: {result['standard_type']}")
        # This line prints out the number in the final result, as requested.
        print(f"Value: {result['relation']} {result['standard_value']} {result['standard_units']}")
        print("--------------------")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an active internet connection and have installed the 'chembl_webresource_client' library (`pip install chembl_webresource_client`).")

if __name__ == '__main__':
    find_binding_affinity()