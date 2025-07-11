# First, ensure you have the necessary library installed:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_binding_affinity():
    """
    Finds and evaluates the binding affinity of a compound to a target protein
    by querying the ChEMBL database.
    """
    # The provided chemical is a known CDK7 inhibitor, SY-1365.
    compound_synonym = 'SY-1365'
    target_name = 'CDK7'
    organism = 'Homo sapiens'

    print(f"Searching for compound synonym: '{compound_synonym}'")
    print(f"Searching for target: '{target_name}' in organism '{organism}'")
    print("-" * 30)

    try:
        # Step 1: Find the ChEMBL ID for the compound 'SY-1365'
        molecule = new_client.molecule
        mol_results = molecule.search(compound_synonym)
        if not mol_results:
            print(f"Error: Compound '{compound_synonym}' not found in ChEMBL.")
            return
        compound_chembl_id = mol_results[0]['molecule_chembl_id']
        print(f"Found Compound ChEMBL ID: {compound_chembl_id}")

        # Step 2: Find the ChEMBL ID for the human CDK7 target
        target = new_client.target
        target_results = target.filter(pref_name__iexact=target_name, organism=organism)
        if not target_results:
            print(f"Error: Target '{target_name}' for '{organism}' not found in ChEMBL.")
            return
        target_chembl_id = target_results[0]['target_chembl_id']
        print(f"Found Target ChEMBL ID: {target_chembl_id}")

        # Step 3: Retrieve bioactivity data (IC50) for the compound-target pair
        activity = new_client.activity
        activity_results = activity.filter(
            molecule_chembl_id=compound_chembl_id,
            target_chembl_id=target_chembl_id,
            standard_type="IC50"  # IC50 is a standard measure of inhibitory concentration
        )

        if not activity_results:
            print(f"\nNo IC50 activity data found for {compound_synonym} against {target_name}.")
            return

        # Step 4: Extract and display the most relevant result
        # We will use the first result, as it is often from a primary publication.
        result = activity_results[0]
        affinity_value = float(result['standard_value'])
        affinity_units = result['standard_units']

        print("\n--- Binding Affinity Data Found ---")
        print(f"Binding Affinity (IC50) = {affinity_value} {affinity_units}")
        print("-" * 30)
        
        # Step 5: Compare with the given answer choices
        print("\nComparing value with the provided ranges:")
        print(f"A. < 0.1 nM")
        print(f"B. 0.1 - 100 nM")
        print(f"C. 0.1 - 100 uM")
        print(f"D. 0.1 - 100 mM")
        print(f"E. > 100 mM")

        if affinity_units == 'nM':
            if 0.1 <= affinity_value <= 100:
                print(f"\nThe value {affinity_value} {affinity_units} falls into range B.")
            elif affinity_value < 0.1:
                print(f"\nThe value {affinity_value} {affinity_units} falls into range A.")
            # Convert to uM range (0.1uM = 100nM, 100uM = 100,000nM)
            elif 100 < affinity_value <= 100000:
                 print(f"\nThe value {affinity_value} {affinity_units} falls into range C.")
            else:
                 print(f"\nThe value {affinity_value} {affinity_units} is outside of the A-C ranges.")
        else:
            print(f"Cannot perform comparison because units are {affinity_units}, not nM.")

    except Exception as e:
        print(f"\nAn error occurred during the search: {e}")
        print("Please ensure you have an active internet connection and the 'chembl_webresource_client' library is installed.")

if __name__ == '__main__':
    find_binding_affinity()
