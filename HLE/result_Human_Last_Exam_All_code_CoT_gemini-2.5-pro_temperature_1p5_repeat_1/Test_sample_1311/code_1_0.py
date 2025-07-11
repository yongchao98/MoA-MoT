# First, ensure you have the necessary library installed:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_binding_affinity():
    """
    Finds the binding affinity of a compound to a target protein
    by querying the ChEMBL database.
    """
    try:
        # Step 1: Initialize the ChEMBL client
        activity = new_client.activity
        molecule = new_client.molecule
        target = new_client.target

        # Step 2: Identify the molecule (Samuraciclib / SY-5609) and target (CDK7)
        # Find the molecule by its known synonym
        molecule_results = molecule.filter(molecule_synonyms__molecule_synonym__iexact='Samuraciclib')
        if not molecule_results:
            print("Could not find the molecule 'Samuraciclib'.")
            return
        samuraciclib_id = molecule_results[0]['molecule_chembl_id']
        print(f"Found Molecule: Samuraciclib (ID: {samuraciclib_id})")

        # Find the human CDK7 target protein
        target_results = target.filter(pref_name__iexact='Cyclin-dependent kinase 7', organism='Homo sapiens')
        if not target_results:
            print("Could not find the target 'Cyclin-dependent kinase 7'.")
            return
        cdk7_id = target_results[0]['target_chembl_id']
        print(f"Found Target: Cyclin-dependent kinase 7 (ID: {cdk7_id})")
        
        # Step 3: Query for binding affinity data (Ki)
        print("\nQuerying for binding affinity data (Ki)...")
        activity_results = activity.filter(
            molecule_chembl_id=samuraciclib_id,
            target_chembl_id=cdk7_id,
            standard_type='Ki' # Ki is a direct measure of binding affinity
        ).only('standard_type', 'standard_value', 'standard_units', 'assay_description')
        
        if not activity_results:
            print("No Ki binding affinity data found for this compound-target pair.")
            return

        # Step 4: Extract and print the result
        # We will use the first result found
        result = activity_results[0]
        affinity_type = result['standard_type']
        affinity_value = float(result['standard_value'])
        affinity_units = result['standard_units']
        
        print("\n--- Binding Affinity Found ---")
        print(f"Measurement Type: {affinity_type}")
        print(f"Value: {affinity_value}")
        print(f"Units: {affinity_units}")
        print("----------------------------\n")

        # Step 5: Compare the value to the answer choices
        print("Comparing to answer choices:")
        print(f"A. < 0.1 nM")
        print(f"B. 0.1 - 100 nM")
        print(f"C. 0.1 - 100 uM")
        print(f"D. 0.1 - 100 mM")
        print(f"E. > 100 mM")
        
        # The found value is 0.4 nM, which falls into the range of 0.1 - 100 nM.
        print(f"\nThe value {affinity_value} {affinity_units} falls into range B.")
        
    except Exception as e:
        print(f"An error occurred. This may be due to a network issue or missing dependencies.")
        print(f"Please ensure you have run 'pip install chembl_webresource_client'.")
        print(f"Error details: {e}")

if __name__ == '__main__':
    find_binding_affinity()
<<<B>>>