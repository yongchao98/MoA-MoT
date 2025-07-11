# First, ensure you have the ChEMBL web resource client installed:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_binding_affinity():
    """
    Finds the binding affinity of a specific compound to CDK7 by querying
    the ChEMBL database.
    """
    try:
        # Initialize clients for different ChEMBL entities
        activity_client = new_client.activity
        molecule_client = new_client.molecule
        target_client = new_client.target

        # Step 1: Find the target protein (CDK7) for Homo sapiens
        print("Searching for Target: Cyclin-dependent kinase 7 (CDK7)...")
        target_search_results = target_client.search('CDK7')
        human_cdk7_target = None
        for result in target_search_results:
            if result['organism'] == 'Homo sapiens':
                human_cdk7_target = result
                break

        if not human_cdk7_target:
            print("Could not find target CDK7 for Homo sapiens.")
            return

        target_chembl_id = human_cdk7_target['target_chembl_id']
        print(f"Found Target: {human_cdk7_target['pref_name']} with ChEMBL ID: {target_chembl_id}\n")

        # Step 2: Find the compound using its InChIKey.
        # The InChIKey for the compound is LJRHSVLXACXBSH-QMMMGPOBSA-N.
        compound_name = "7-dimethylphosphoryl-3-[2-[[(3~{S})-6,6-dimethylpiperidin-3-yl]amino]-5-(trifluoromethyl)pyrimidin-4-yl]-1~{H}-indole-6-carbonitrile"
        compound_inchi_key = 'LJRHSVLXACXBSH-QMMMGPOBSA-N'
        print(f"Searching for Compound by its unique identifier (InChIKey)...")

        molecule_search_results = molecule_client.filter(inchi_key=compound_inchi_key)
        if not molecule_search_results:
            print(f"Could not find compound with InChIKey: {compound_inchi_key}")
            return

        compound_info = molecule_search_results[0]
        molecule_chembl_id = compound_info['molecule_chembl_id']
        print(f"Found Compound (prefers name: {compound_info.get('pref_name', 'N/A')}) with ChEMBL ID: {molecule_chembl_id}\n")

        # Step 3: Retrieve binding affinity data (Ki, IC50)
        print("Querying for binding affinity data...\n")
        activities = activity_client.filter(
            target_chembl_id=target_chembl_id,
            molecule_chembl_id=molecule_chembl_id,
            standard_type__in=['Ki', 'IC50'] # Look for Ki or IC50 values
        ).only(['type', 'standard_units', 'standard_value', 'relation'])

        if not activities:
            print("No Ki or IC50 binding affinity data found for this compound-target pair in ChEMBL.")
            return

        print("Found the following binding affinity values:")
        found_in_range = False
        for act in activities:
            value = act['standard_value']
            units = act['standard_units']
            relation = act['relation']
            # We are interested in nanomolar (nM) range values
            if value and units == 'nM':
                print(f"-> Activity Type: {act['type']}, Value: {relation}{value} {units}")
                # Check if the value falls in the B range (0.1 - 100 nM)
                if float(value) >= 0.1 and float(value) <= 100:
                    found_in_range = True

        if found_in_range:
             print("\nConclusion: The reported binding affinities fall within the 0.1 - 100 nM range.")
        else:
             print("\nConclusion: While values were found, they do not fall within the 0.1 - 100 nM range.")


    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an internet connection and the 'chembl_webresource_client' library is installed ('pip install chembl_webresource_client').")

# Execute the function
if __name__ == "__main__":
    find_binding_affinity()