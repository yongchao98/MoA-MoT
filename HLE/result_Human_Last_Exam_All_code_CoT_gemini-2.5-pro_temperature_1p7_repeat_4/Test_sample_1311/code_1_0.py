# The user will need to install the chembl_webresource_client library
# You can do this by running: pip install chembl_webresource_client pandas

from chembl_webresource_client.new_client import new_client
import pandas as pd
import sys

def find_binding_affinity():
    """
    Finds the binding affinity of Samuraciclib to human CDK7 using the ChEMBL database.
    """
    try:
        # Initialize the ChEMBL new_client
        activity = new_client.activity
        target = new_client.target
        molecule = new_client.molecule
    except Exception as e:
        print(f"Failed to connect to ChEMBL API. Please check your internet connection.", file=sys.stderr)
        print("You may also need to install the client library: pip install chembl_webresource_client pandas", file=sys.stderr)
        return

    # --- Step 1: Find the target protein (CDK7) ---
    print("Step 1: Searching for the target protein 'Cyclin-dependent kinase 7'...")
    # Search for the target protein by its preferred name
    target_search = target.search('Cyclin-dependent kinase 7')
    if not target_search:
        print("Error: Target 'Cyclin-dependent kinase 7' not found.", file=sys.stderr)
        return
        
    targets_df = pd.DataFrame.from_records(target_search)
    # Filter for the human protein
    human_cdk7_target = targets_df[(targets_df['organism'] == 'Homo sapiens') & (targets_df['target_type'] == 'SINGLE PROTEIN')]
    
    if human_cdk7_target.empty:
        print("Error: Human 'Cyclin-dependent kinase 7' not found.", file=sys.stderr)
        return
        
    target_chembl_id = human_cdk7_target.iloc[0]['target_chembl_id']
    print(f"Found target: {human_cdk7_target.iloc[0]['pref_name']} (Homo sapiens), ChEMBL ID: {target_chembl_id}\n")

    # --- Step 2: Find the molecule (Samuraciclib) ---
    compound_name = "Samuraciclib"
    print(f"Step 2: Searching for the molecule '{compound_name}'...")
    # The provided IUPAC name is: 7-dimethylphosphoryl-3-[2-[[(3~{S})-6,6-dimethylpiperidin-3-yl]amino]-5-(trifluoromethyl)pyrimidin-4-yl]-1~{H}-indole-6-carbonitrile
    molecule_search = molecule.search(compound_name)
    if not molecule_search:
        print(f"Error: Molecule '{compound_name}' not found.", file=sys.stderr)
        return

    molecules_df = pd.DataFrame.from_records(molecule_search)
    molecule_chembl_id = molecules_df.iloc[0]['molecule_chembl_id']
    print(f"Found molecule: {molecules_df.iloc[0]['pref_name']}, ChEMBL ID: {molecule_chembl_id}\n")

    # --- Step 3: Retrieve binding affinity data ---
    print(f"Step 3: Fetching binding affinity data for {compound_name} and {human_cdk7_target.iloc[0]['pref_name']}...")
    # Filter for activity data linking the molecule and the target
    # We look for standard binding affinity measures like Ki, Kd, IC50
    activity_data = activity.filter(
        target_chembl_id=target_chembl_id,
        molecule_chembl_id=molecule_chembl_id,
        standard_type__in=['Ki', 'IC50', 'Kd']
    ).only(['standard_type', 'standard_relation', 'standard_value', 'standard_units', 'assay_description'])

    if not activity_data:
        print("No binding affinity data (Ki, IC50, Kd) found for this pair in the ChEMBL database.", file=sys.stderr)
        return
    
    print("Found the following relevant bioactivity records:")
    
    affinity_values_nm = []
    for record in activity_data:
        # We only consider data with units of 'nM' as it's most common for potent inhibitors
        if record['standard_units'] == 'nM' and record['standard_value'] is not None:
            value = float(record['standard_value'])
            affinity_values_nm.append(value)
            print(f"- Assay: '{record['assay_description']}'")
            # This fulfills the "output each number in the final equation" requirement
            print(f"  Result: {record['standard_type']} {record['standard_relation']} {record['standard_value']} {record['standard_units']}")

    if not affinity_values_nm:
        print("\nCould not find any binding affinity values measured in nM.", file=sys.stderr)
        return

    # --- Step 4: Analyze the results ---
    print("\nStep 4: Analyzing the results...")
    # Check if all found values fall within a specific range
    is_in_range_B = all(0.1 <= val <= 100 for val in affinity_values_nm)
    
    print("All retrieved affinity values are in the nanomolar (nM) range.")
    print("The answer choices are:")
    print("A. < 0.1 nM")
    print("B. 0.1 - 100 nM")
    print("C. 0.1 - 100 uM")
    print("D. 0.1 - 100 mM")
    print("E. > 100 mM")
    
    if is_in_range_B:
        print("\nConclusion: All recorded binding affinities fall within the 0.1 - 100 nM range.")
    else:
        # This case is unlikely based on prior knowledge, but included for robustness
        print("\nConclusion: The recorded binding affinities fall into different ranges.")

if __name__ == "__main__":
    find_binding_affinity()
<<<B>>>