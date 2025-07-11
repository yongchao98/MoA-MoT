# First, you may need to install the ChEMBL web resource client library
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_interacting_molecules(target_chembl_id):
    """
    Finds and prints small molecules that interact with a given ChEMBL target ID.
    The search is filtered for potent binders (pChEMBL value >= 5.0).
    """
    print(f"Searching for small molecule ligands for target: {target_chembl_id}\n")
    
    try:
        # Initialize the client for accessing bioactivities
        activity = new_client.activity

        # Query for activities related to the target ID
        # Filter for pChEMBL value >= 5 to get reasonably potent compounds
        # (pChEMBL of 5 corresponds to an IC50/EC50/Ki of 10,000 nM or 10 uM)
        res = activity.filter(target_chembl_id=target_chembl_id, pchembl_value__gte=5).only(['molecule_chembl_id', 'pchembl_value', 'standard_type'])

        if not res:
            print(f"No interacting molecules found for {target_chembl_id} with a pChEMBL value >= 5.0.")
            return

        # Use a set to store unique molecule IDs to avoid duplicates
        interacting_molecule_ids = set()
        for item in res:
            if item['molecule_chembl_id']:
                interacting_molecule_ids.add(item['molecule_chembl_id'])

        print(f"Found {len(interacting_molecule_ids)} unique small molecules interacting with {target_chembl_id}:")
        
        # Print each molecule ID
        for molecule_id in sorted(list(interacting_molecule_ids)):
            print(molecule_id)

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure the 'chembl_webresource_client' library is installed (`pip install chembl_webresource_client`) and you have an internet connection.")


if __name__ == '__main__':
    target_id = "CHEMBL4803817"
    find_interacting_molecules(target_id)