# First, ensure you have the ChEMBL client library installed:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_binding_affinity():
    """
    Finds the binding affinity (IC50) of a molecule to a target protein
    by querying the ChEMBL database.
    """
    try:
        # Initialize clients for different ChEMBL resources
        molecule = new_client.molecule
        target = new_client.target
        activity = new_client.activity

        # --- Step 1: Identify the molecule ---
        # The molecule is SY-1365
        molecule_name = 'SY-1365'
        print(f"Searching for molecule: {molecule_name}")
        mol_results = molecule.filter(synonyms__icontains=molecule_name)
        if not mol_results:
            print(f"Could not find molecule '{molecule_name}' in ChEMBL.")
            return

        target_molecule = mol_results[0]
        molecule_chembl_id = target_molecule['molecule_chembl_id']
        print(f"Found molecule ChEMBL ID: {molecule_chembl_id}\n")

        # --- Step 2: Identify the target protein ---
        # The target is human Cyclin-dependent kinase 7 (CDK7)
        target_name = 'Cyclin-dependent kinase 7'
        target_organism = 'Homo sapiens'
        print(f"Searching for target: {target_name} ({target_organism})")
        target_results = target.filter(
            pref_name__iexact=target_name,
            target_organism__iexact=target_organism
        )
        if not target_results:
            print(f"Could not find target '{target_name}' for organism '{target_organism}' in ChEMBL.")
            return
            
        target_protein = target_results[0]
        target_chembl_id = target_protein['target_chembl_id']
        print(f"Found target ChEMBL ID: {target_chembl_id}\n")

        # --- Step 3: Retrieve bioactivity data ---
        # We are looking for IC50 values
        print(f"Querying for IC50 values for {molecule_chembl_id} on {target_chembl_id}...")
        activity_results = activity.filter(
            molecule_chembl_id=molecule_chembl_id,
            target_chembl_id=target_chembl_id,
            standard_type="IC50"
        )

        if not activity_results:
            print("No IC50 activity data found for this molecule-target pair.")
            return

        # --- Step 4: Print the results and determine the range ---
        print("\nFound binding affinity data:")
        for res in activity_results:
            value = res['standard_value']
            units = res['standard_units']
            if value and units == 'nM':
                value_float = float(value)
                print(f"IC50 = {value_float} {units}")

                # Compare with the answer choices
                if 0.1 <= value_float <= 100:
                    print(f"This value falls within the range B: 0.1 - 100 nM")
                elif value_float < 0.1:
                    print(f"This value falls within the range A: < 0.1 nM")
                else:
                    print("This value does not fall into the primary nM ranges listed.")
                # We only need one representative value to answer the question
                return

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an internet connection and the 'chembl_webresource_client' library is installed.")

if __name__ == '__main__':
    find_binding_affinity()
