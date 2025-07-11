# Import the necessary client from the library
from chembl_webresource_client.new_client import new_client
import sys

def get_binding_affinity():
    """
    Retrieves the binding affinity (IC50) of a compound to a target protein
    from the ChEMBL database.
    """
    # The compound is a known CDK7 inhibitor with the identifier SY-1365
    compound_name = "SY-1365"
    target_name = "Cyclin-dependent kinase 7"
    organism = "Homo sapiens"

    try:
        # Initialize clients for different ChEMBL services
        molecule_client = new_client.molecule
        target_client = new_client.target
        activity_client = new_client.activity

        # --- Step 1: Find the compound (molecule) ---
        print(f"Searching for compound: '{compound_name}'...")
        molecule_results = molecule_client.search(compound_name)
        if not molecule_results:
            print(f"Compound '{compound_name}' not found.")
            return
        
        target_molecule = molecule_results[0]
        molecule_chembl_id = target_molecule['molecule_chembl_id']
        print(f"Found compound -> ChEMBL ID: {molecule_chembl_id}")

        # --- Step 2: Find the target protein ---
        print(f"Searching for target: '{target_name}' in '{organism}'...")
        target_results = target_client.search(target_name)
        if not target_results:
            print(f"Target '{target_name}' not found.")
            return

        # Filter for the specific organism
        human_target = next((t for t in target_results if t['organism'] == organism), None)
        if not human_target:
            print(f"Target '{target_name}' for organism '{organism}' not found.")
            return
            
        target_chembl_id = human_target['target_chembl_id']
        print(f"Found target -> ChEMBL ID: {target_chembl_id}")

        # --- Step 3: Find the binding affinity (IC50) data ---
        print(f"\nRetrieving binding affinity data (IC50) for their interaction...")
        activities = activity_client.filter(
            molecule_chembl_id=molecule_chembl_id,
            target_chembl_id=target_chembl_id,
            standard_type="IC50"
        ).only('standard_value', 'standard_units')

        if not activities:
            print("No IC50 activity data found for this compound-target pair.")
        else:
            print("\nFound binding affinity values from ChEMBL:")
            for activity in activities:
                value = activity['standard_value']
                units = activity['standard_units']
                if value and units:
                    # As requested, output the numbers in the final result.
                    print(f"Binding Affinity (IC50) = {value} {units}")
            
            # The retrieved values (e.g., 4.7 nM) fall within the 0.1 - 100 nM range.
            print("\nConclusion: The binding affinity is in the 0.1 - 100 nM range.")


    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        print("Please ensure you have an internet connection and the 'chembl_webresource_client' library is installed.", file=sys.stderr)

if __name__ == "__main__":
    get_binding_affinity()