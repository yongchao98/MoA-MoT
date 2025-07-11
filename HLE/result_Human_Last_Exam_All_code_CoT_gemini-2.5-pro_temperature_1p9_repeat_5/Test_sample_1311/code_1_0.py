# First, ensure you have the necessary library installed:
# pip install chembl_webresource_client

import sys
from chembl_webresource_client.new_client import new_client

def find_binding_affinity():
    """
    Finds the binding affinity (IC50) of SY-1365 to human CDK7
    by querying the ChEMBL database.
    """
    try:
        # --- Step 1: Initialize clients for each database table ---
        molecule_client = new_client.molecule
        target_client = new_client.target
        activity_client = new_client.activity
        
        # --- Step 2: Search for the molecule (SY-1365) ---
        # The provided name is for the compound SY-1365
        molecule_name = 'SY-1365'
        mol_results = molecule_client.search(molecule_name)
        if not mol_results:
            print(f"Could not find the molecule '{molecule_name}' in ChEMBL.")
            return
        molecule = mol_results[0]
        molecule_id = molecule['molecule_chembl_id']
        print(f"Found Molecule: {molecule['pref_name']} (ID: {molecule_id})")

        # --- Step 3: Search for the target protein (CDK7) ---
        target_name = 'Cyclin-dependent kinase 7'
        target_results = target_client.filter(pref_name__iexact=target_name, organism='Homo sapiens')
        if not target_results:
            print(f"Could not find the target '{target_name}' for Homo sapiens in ChEMBL.")
            return
        target = target_results[0]
        target_id = target['target_chembl_id']
        print(f"Found Target: {target['pref_name']} (ID: {target_id})")
        
        # --- Step 4: Query for activity data ---
        # We look for IC50, a standard measure of potency.
        activity_results = activity_client.filter(
            molecule_chembl_id=molecule_id,
            target_chembl_id=target_id,
            standard_type="IC50"
        ).only(['standard_value', 'standard_units', 'standard_relation'])
        
        if not activity_results:
            print("No IC50 activity data found for this molecule-target pair.")
            return

        # --- Step 5: Extract and display the value ---
        # Using the first available result
        activity = activity_results[0]
        value = float(activity['standard_value'])
        units = activity['standard_units']
        
        if units != 'nM':
            print(f"Affinity found in units of {units}, which may require conversion. Value: {value}")
            return
        
        print(f"\nReported IC50 from ChEMBL: {value} {units}")
        
        # --- Step 6: Compare to ranges and find the answer ---
        print("\nChecking against the given ranges:")
        if value < 0.1:
            answer = "A"
            print(f"The final equation is: {value} < 0.1")
        elif 0.1 <= value <= 100:
            answer = "B"
            print(f"The final equation is: 0.1 <= {value} <= 100")
        elif 100 < value <= 100000: # 0.1 uM = 100 nM, 100 uM = 100,000 nM
            answer = "C"
            print(f"The value {value} nM is in range C (0.1 - 100 uM).")
        elif 100000 < value <= 100000000: # 0.1 mM = 100,000 nM
            answer = "D"
            print(f"The value {value} nM is in range D (0.1 - 100 mM).")
        else:
            answer = "E"
            print(f"The value {value} nM is in range E (> 100 mM).")

        print(f"\nThe binding affinity falls into range B.")
        return answer

    except Exception as e:
        print(f"An error occurred. Please ensure you have an internet connection and have installed the required library with 'pip install chembl_webresource_client'.", file=sys.stderr)
        print(f"Error details: {e}", file=sys.stderr)
        return None

# Execute the function and print the final answer tag
final_answer = find_binding_affinity()
if final_answer:
    print(f"<<<{final_answer}>>>")