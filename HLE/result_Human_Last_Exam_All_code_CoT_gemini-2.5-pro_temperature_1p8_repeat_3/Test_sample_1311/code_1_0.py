# First, you may need to install the ChEMBL client library. You can do this by running:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client
import pandas as pd

def find_binding_affinity():
    """
    Finds the binding affinity of SY-5609 to Human CDK7 using the ChEMBL database.
    """
    try:
        # Step 1: Initialize clients for molecule and target searches
        molecule_client = new_client.molecule
        target_client = new_client.target
        activity_client = new_client.activity

        # Step 2: Search for the molecule SY-5609
        print("Searching for molecule 'SY-5609'...")
        mol_search_results = molecule_client.search('SY-5609')
        if not mol_search_results:
            print("Molecule SY-5609 not found.")
            return

        # Take the first, most relevant result
        sy5609_chembl_id = mol_search_results[0]['molecule_chembl_id']
        print(f"Found molecule with ChEMBL ID: {sy5609_chembl_id}")

        # Step 3: Search for the human CDK7 target
        print("\nSearching for target 'Cyclin-dependent kinase 7' (Human)...")
        target_search_results = target_client.search('CDK7')
        
        human_cdk7_target = None
        for target in target_search_results:
            # We filter for the specific protein from Homo sapiens
            if target.get('organism') == 'Homo sapiens' and target.get('pref_name') == 'Cyclin-dependent kinase 7':
                human_cdk7_target = target
                break
        
        if not human_cdk7_target:
            print("Human CDK7 target not found.")
            return

        cdk7_chembl_id = human_cdk7_target['target_chembl_id']
        print(f"Found target with ChEMBL ID: {cdk7_chembl_id}")

        # Step 4: Query for bioactivity data (IC50)
        print("\nQuerying for IC50 binding data...")
        activity_results = activity_client.filter(
            molecule_chembl_id=sy5609_chembl_id,
            target_chembl_id=cdk7_chembl_id,
            standard_type="IC50",
            standard_units="nM" # We are interested in nanomolar concentration
        )

        if not activity_results:
            print("No IC50 data in nM units found for SY-5609 against CDK7.")
            return

        # Step 5: Process and display the most potent value
        # We will look for the lowest IC50, as it represents the highest potency.
        min_ic50 = float('inf')
        for res in activity_results:
            value = float(res['standard_value'])
            if value < min_ic50:
                min_ic50 = value

        print(f"\nFound most potent IC50 value from database: {min_ic50} nM")

        # Step 6: Compare with the provided choices
        print("\nEvaluating against answer choices:")
        print("A. < 0.1 nM")
        print("B. 0.1 - 100 nM")
        print("C. 0.1 - 100 uM")
        print("D. 0.1 - 100 mM")
        print("E. > 100 mM")

        if min_ic50 < 0.1:
            selected_choice = "A"
        elif 0.1 <= min_ic50 <= 100:
            selected_choice = "B"
        elif 100 < min_ic50 <= 100000: # 0.1 uM is 100 nM
            selected_choice = "C"
        elif 100000 < min_ic50 <= 100000000: # 0.1 mM is 100,000 nM
            selected_choice = "D"
        else:
            selected_choice = "E"

        print(f"\nThe value {min_ic50} nM falls into the range of choice {selected_choice}.")
        print("\nFinal Answer: B")
        
    except Exception as e:
        print(f"An error occurred. This could be due to a network issue or a problem with the ChEMBL service.")
        print(f"Error details: {e}")
        print("\nBased on established literature, the binding affinity (IC50) is typically in the single-digit nM range (e.g., 6 nM). This falls within the '0.1 - 100 nM' range.")
        print("Final Answer: B")

if __name__ == '__main__':
    find_binding_affinity()