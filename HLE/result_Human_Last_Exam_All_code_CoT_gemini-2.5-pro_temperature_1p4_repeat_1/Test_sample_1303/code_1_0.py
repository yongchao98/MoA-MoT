# First, ensure you have the necessary library installed:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_interacting_molecules(target_id):
    """
    Finds and prints small molecules that interact with a given ChEMBL target ID.
    """
    try:
        # Initialize the ChEMBL web resource clients
        target = new_client.target
        activity = new_client.activity

        # --- Step 1: Fetch target information ---
        print(f"Searching for target with ChEMBL ID: {target_id}")
        target_info = target.get(target_id)
        
        if not target_info:
            print(f"Could not find a target with the ID {target_id}.")
            return

        target_name = target_info.get('pref_name', 'N/A')
        print(f"Found Target: {target_name} ({target_id})\n")

        # --- Step 2: Fetch associated activities ---
        print(f"Finding small molecules that interact with {target_name}...")
        # We filter for activities associated with the target and select only the molecule ID field
        activities = activity.filter(target_chembl_id=target_id).only(['molecule_chembl_id'])

        if not activities:
            print(f"No interacting small molecules found for {target_name} in the ChEMBL database.")
            return

        # --- Step 3: Extract and print unique molecule IDs ---
        # Use a set to store unique molecule IDs to avoid duplicates
        interacting_molecule_ids = set()
        for act in activities:
            if 'molecule_chembl_id' in act and act['molecule_chembl_id']:
                interacting_molecule_ids.add(act['molecule_chembl_id'])
        
        if interacting_molecule_ids:
            print("\nList of interacting small molecule ChEMBL IDs:")
            # Sort the IDs for consistent output
            sorted_ids = sorted(list(interacting_molecule_ids))
            for molecule_id in sorted_ids:
                print(molecule_id)
        else:
             print("Found activity entries, but could not extract molecule ChEMBL IDs.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an active internet connection and have installed the 'chembl_webresource_client' library.")
        print("You can install it using: pip install chembl_webresource_client")

if __name__ == "__main__":
    # The ChEMBL ID for the target of interest
    target_chembl_id = 'CHEMBL4803817'
    find_interacting_molecules(target_chembl_id)
