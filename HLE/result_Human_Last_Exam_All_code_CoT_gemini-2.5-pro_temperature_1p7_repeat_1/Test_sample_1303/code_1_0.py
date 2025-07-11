try:
    # First, ensure you have the library installed:
    # pip install chembl_webresource_client
    from chembl_webresource_client.new_client import new_client

    # The ChEMBL ID for the target protein (SARS-CoV-2 protein E)
    target_chembl_id = "CHEMBL4803817"

    # Access the activity data from ChEMBL
    activity = new_client.activity

    # Search for activities related to the target ID.
    # We filter for activities with a pChEMBL value to ensure a measured interaction.
    print(f"Searching for small molecules that interact with target: {target_chembl_id}...")
    res = activity.filter(target_chembl_id=target_chembl_id, pchembl_value__isnull=False)

    if res:
        # Use a set to store unique molecule IDs to avoid duplicates
        molecule_ids = set()
        for item in res:
            molecule_ids.add(item['molecule_chembl_id'])

        print(f"\nFound {len(molecule_ids)} unique small molecules interacting with {target_chembl_id}:")
        # Print each molecule's ChEMBL ID
        for mol_id in sorted(list(molecule_ids)):
            print(mol_id)
    else:
        print(f"No interacting small molecules with a pChEMBL value were found for {target_chembl_id}.")

except ImportError:
    print("The 'chembl_webresource_client' library is not installed.")
    print("Please install it using: pip install chembl_webresource_client")
except Exception as e:
    print(f"An error occurred: {e}")
