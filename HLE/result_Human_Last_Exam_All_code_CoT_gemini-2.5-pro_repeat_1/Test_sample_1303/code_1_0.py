# First, you may need to install the ChEMBL web resource client
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_interacting_molecules(target_chembl_id):
    """
    Finds and prints small molecules that interact with a given ChEMBL target ID.
    """
    try:
        # Get the target resource
        target = new_client.target
        # Get the activity resource
        activity = new_client.activity

        print(f"Searching for small molecule interactions with target: {target_chembl_id}\n")

        # Query for activities related to the target
        activities = activity.filter(target_chembl_id=target_chembl_id).only(
            ['molecule_chembl_id', 'molecule_type', 'standard_type', 'standard_value', 'standard_units']
        )

        interacting_molecules = set()
        for act in activities:
            # Filter for small molecules with a standard activity value
            if act['molecule_type'] == 'Small molecule' and act['standard_value'] is not None:
                interacting_molecules.add(act['molecule_chembl_id'])

        if interacting_molecules:
            print(f"Found {len(interacting_molecules)} unique small molecules interacting with {target_chembl_id}:")
            for molecule_id in sorted(list(interacting_molecules)):
                print(molecule_id)
        else:
            print(f"No small molecule interactions with a standard value were found for {target_chembl_id}.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure the ChEMBL ID is correct and you have an internet connection.")

if __name__ == "__main__":
    target_id = "CHEMBL4803817"
    find_interacting_molecules(target_id)