# First, ensure you have the ChEMBL web resource client installed.
# You can install it using pip:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client
import sys

def find_interacting_molecules(component_id):
    """
    Finds small molecules that interact with a target specified by a ChEMBL Component ID.
    """
    try:
        # Get the target information from the component ID
        target_component = new_client.target_component
        target_info = target_component.filter(component_chembl_id=component_id).only('target_chembl_id')

        if not target_info:
            print(f"Error: Could not find a target for component ID '{component_id}'.")
            return None

        target_id = target_info[0]['target_chembl_id']
        print(f"Found Target ChEMBL ID: {target_id}")
        print("-" * 30)

        # Query for activities involving this target and small molecules
        activity = new_client.activity
        activities_query = activity.filter(
            target_chembl_id=target_id,
            molecule_type='Small molecule'
        )

        # Get the total count of interacting small molecules
        total_count = activities_query.count()

        if total_count == 0:
            print(f"No small molecule interactions found for target {target_id}.")
            return 0
        
        print(f"Found a total of {total_count} small molecule interactions.")
        print("Here is a sample of the first 20 unique interacting molecules:")
        print("-" * 30)

        # Retrieve a sample of unique molecule IDs
        seen_molecules = set()
        # Limit the initial fetch to avoid pulling too much data for just a sample
        for act in activities_query[:500]: 
            if len(seen_molecules) >= 20:
                break
            if act['molecule_chembl_id']:
                seen_molecules.add(act['molecule_chembl_id'])

        for molecule_id in sorted(list(seen_molecules)):
            print(molecule_id)
            
        return total_count

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        return None

if __name__ == "__main__":
    chembl_id = "CHEMBL4803817"
    print(f"Searching for small molecules interacting with target related to {chembl_id}...\n")
    total_molecules = find_interacting_molecules(chembl_id)
    if total_molecules is not None:
        # This is where the final answer format is used.
        # The number represents the total count of interacting small molecules.
        print(f"\nFinal Answer formatted as requested:")
        print(f"<<<{total_molecules}>>>")
