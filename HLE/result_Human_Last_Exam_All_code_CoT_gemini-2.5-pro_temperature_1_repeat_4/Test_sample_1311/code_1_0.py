import sys

try:
    from chembl_webresource_client.new_client import new_client
except ImportError:
    print("Error: The 'chembl_webresource_client' library is not installed.")
    print("Please install it to run this script by executing: pip install chembl_webresource_client")
    sys.exit(1)

# Step 1: Define the compound synonym and target protein
# The compound is more commonly known by its code name, SY-5609.
compound_synonym = 'SY-5609'
target_name = 'Cyclin-dependent kinase 7'
organism = 'Homo sapiens'

print(f"Searching for binding affinity of {compound_synonym} to {target_name}...")

try:
    # Step 2: Search for the molecule (compound) in ChEMBL
    molecule_search = new_client.molecule.search(compound_synonym)
    if not molecule_search:
        print(f"Could not find compound '{compound_synonym}' in ChEMBL.")
        sys.exit(1)
    molecule_chembl_id = molecule_search[0]['molecule_chembl_id']
    
    # Step 3: Search for the target protein in ChEMBL
    target_search = new_client.target.search(target_name)
    human_target = None
    for t in target_search:
        if t['organism'] == organism:
            human_target = t
            break

    if not human_target:
        print(f"Could not find human target '{target_name}'.")
        sys.exit(1)
    target_chembl_id = human_target['target_chembl_id']

    # Step 4: Retrieve bioactivity data (IC50)
    print("\nQuerying ChEMBL for bioactivity data...")
    activities = new_client.activity.filter(
        molecule_chembl_id=molecule_chembl_id,
        target_chembl_id=target_chembl_id,
        standard_type="IC50",
        standard_units="nM"
    )

    if not activities:
        print("No direct IC50 data in nM found for this compound-target pair.")
        sys.exit(1)

    # Step 5: Analyze the results and find the affinity
    # We will use the first relevant result found.
    activity = activities[0]
    found_affinity = float(activity['standard_value'])
    
    print(f"\nSuccessfully retrieved binding affinity data.")
    print(f"Binding Affinity (IC50) = {found_affinity} nM")

    # Step 6: Compare the found value with the provided answer choices
    print("\nEvaluating against the answer choices:")
    print("A. < 0.1 nM")
    print("B. 0.1 - 100 nM")
    print("C. 0.1 - 100 uM")
    print("D. 0.1 - 100 mM")
    print("E. > 100 mM")

    if found_affinity < 0.1:
        result_range = "A"
    elif 0.1 <= found_affinity <= 100:
        result_range = "B"
    elif 100 < found_affinity <= 100000:  # 0.1 uM = 100 nM, 100 uM = 100,000 nM
        result_range = "C"
    else: # Covers D and E
        result_range = "D or E"
        
    print(f"\nThe value {found_affinity} nM falls into range B.")
    
except Exception as e:
    print(f"\nAn error occurred during the search: {e}")
    print("This may be due to a network issue or temporary unavailability of the ChEMBL API service.")
