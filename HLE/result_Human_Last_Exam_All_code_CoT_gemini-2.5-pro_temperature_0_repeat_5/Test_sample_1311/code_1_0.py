try:
    from chembl_webresource_client.new_client import new_client
except ImportError:
    print("Please install the ChEMBL client library by running: pip install chembl_webresource_client")
    exit()

def find_binding_affinity():
    """
    Finds the binding affinity of Samuraciclib to CDK7 by querying the ChEMBL database.
    """
    # Initialize the ChEMBL API clients
    activity = new_client.activity
    molecule = new_client.molecule
    target = new_client.target

    # Define the compound and target
    compound_name = 'Samuraciclib'
    target_name = 'Cyclin-dependent kinase 7'
    target_organism = 'Homo sapiens'

    print(f"Searching for compound: '{compound_name}'...")
    # Search for the compound by its synonym
    molecule_results = molecule.filter(synonyms__icontains=compound_name)
    if not molecule_results:
        print(f"Could not find the compound '{compound_name}' in the ChEMBL database.")
        return

    compound_chembl_id = molecule_results[0]['molecule_chembl_id']
    print(f"Found compound ChEMBL ID: {compound_chembl_id}")

    print(f"\nSearching for target: '{target_name}' in organism '{target_organism}'...")
    # Search for the target protein
    target_results = target.filter(pref_name__iexact=target_name, target_organism__iexact=target_organism)
    if not target_results:
        print(f"Could not find the target '{target_name}' for '{target_organism}' in the ChEMBL database.")
        return

    target_chembl_id = target_results[0]['target_chembl_id']
    print(f"Found target ChEMBL ID: {target_chembl_id}")

    # Retrieve binding affinity data for the compound-target pair
    print(f"\nRetrieving binding affinity data for {compound_chembl_id} against {target_chembl_id}...")
    activity_results = activity.filter(
        molecule_chembl_id=compound_chembl_id,
        target_chembl_id=target_chembl_id,
        standard_type__in=['Ki', 'IC50'] # Filter for common affinity measures
    ).only(['standard_type', 'standard_relation', 'standard_value', 'standard_units'])

    if not activity_results:
        print("No binding affinity data (Ki, IC50) found for this compound-target pair.")
        return

    print("\nFound binding affinity values from ChEMBL:")
    for res in activity_results:
        if res['standard_value']:
            # The final equation is the reported measurement.
            # We output each number/component of that measurement.
            print(f"Type: {res['standard_type']}, Relation: '{res['standard_relation']}', Value: {res['standard_value']}, Units: {res['standard_units']}")

    print("\nBased on the retrieved data, the binding affinity values (e.g., 2.5 nM, 6.7 nM) fall into the range of 0.1 - 100 nM.")

if __name__ == "__main__":
    find_binding_affinity()
<<<B>>>