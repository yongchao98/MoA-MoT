# First, ensure you have the ChEMBL web resource client installed:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

# The ChEMBL ID for the target protein
target_chembl_id = 'CHEMBL4803817'

try:
    # Get the activity client
    activity = new_client.activity

    # Search for activities related to the target
    # We filter for:
    # - The specific target ChEMBL ID
    # - Molecules of type 'Small molecule'
    # - A pChEMBL value >= 5 to find reasonably potent interactions
    res = activity.filter(
        target_chembl_id=target_chembl_id,
        pchembl_value__gte=5,
        molecule_type='Small molecule'
    ).only(['molecule_chembl_id', 'pchembl_value'])

    if res:
        print(f"Found small molecules interacting with {target_chembl_id} (pChEMBL >= 5):")
        # Use a set to store unique molecule IDs
        unique_molecule_ids = set()
        for r in res:
            unique_molecule_ids.add(r['molecule_chembl_id'])
        
        # Print each unique molecule ID
        for molecule_id in sorted(list(unique_molecule_ids)):
            print(molecule_id)
        
        print(f"\n<<<Found {len(unique_molecule_ids)} unique small molecules.>>>")

    else:
        print(f"No small molecule interactions with a pChEMBL value >= 5 found for {target_chembl_id}.")
        print("<<<No interacting small molecules found.>>>")

except Exception as e:
    print(f"An error occurred: {e}")
    print("Please ensure you have an internet connection and the 'chembl_webresource_client' library is installed ('pip install chembl_webresource_client').")
    print("<<<Error occurred during execution.>>>")
