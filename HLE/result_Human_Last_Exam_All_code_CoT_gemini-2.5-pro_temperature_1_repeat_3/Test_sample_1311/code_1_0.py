# First, ensure you have the necessary library installed.
# You can install it by running: pip install chembl_webresource_client
import sys
try:
    from chembl_webresource_client.new_client import new_client
except ImportError:
    print("The 'chembl_webresource_client' library is not installed.")
    print("Please install it by running: pip install chembl_webresource_client")
    sys.exit(1)

def find_binding_affinity():
    """
    Finds the binding affinity (IC50) of a compound to a target protein
    by querying the ChEMBL database.
    """
    # Define search terms based on the query.
    # The compound is 7-dimethylphosphoryl-3-[2-[[(3~{S})-6,6-dimethylpiperidin-3-yl]amino]-5-(trifluoromethyl)pyrimidin-4-yl]-1~{H}-indole-6-carbonitrile
    # This compound is also known as SY-5609.
    # We will use its standard InChIKey for a precise search.
    compound_inchi_key = 'YZWMYRIQJPVLRG-CYBMUJFWSA-N'
    target_name = 'Cyclin-dependent kinase 7'
    organism = 'Homo sapiens'

    # Initialize the APIs
    molecule_api = new_client.molecule
    target_api = new_client.target
    activity_api = new_client.activity

    try:
        # Step 1: Find the molecule's ChEMBL ID using its InChIKey
        mol_results = molecule_api.filter(inchi_key=compound_inchi_key).only('molecule_chembl_id', 'pref_name')
        if not mol_results:
            print(f"Error: Could not find compound with InChIKey: {compound_inchi_key}")
            return

        molecule = mol_results[0]
        mol_chembl_id = molecule['molecule_chembl_id']
        mol_name = molecule.get('pref_name', 'Unknown Compound Name')

        # Step 2: Find the target's ChEMBL ID
        target_results = target_api.search(target_name).filter(organism=organism).only('target_chembl_id', 'pref_name')
        if not target_results:
            print(f"Error: Could not find target '{target_name}' for organism '{organism}'")
            return

        target = target_results[0]
        target_chembl_id = target['target_chembl_id']
        target_pref_name = target['pref_name']

        print(f"Querying ChEMBL database for binding affinity...")
        print(f"Compound: {mol_name} ({mol_chembl_id})")
        print(f"Target: {target_pref_name} ({target_chembl_id})")
        print("-" * 40)

        # Step 3: Query for activity data (IC50)
        activities = activity_api.filter(
            molecule_chembl_id=mol_chembl_id,
            target_chembl_id=target_chembl_id,
            standard_type="IC50",
            relation='=' # We want an exact value
        ).only('standard_value', 'standard_units')

        if not activities:
            print(f"No exact IC50 activity data found for {mol_name} against {target_pref_name}.")
            return

        # Step 4: Process and print the result
        # We prioritize results in nM units as they are common for affinity.
        found_nM_value = False
        for act in activities:
            if act['standard_units'] == 'nM':
                ic50_value = float(act['standard_value'])
                print(f"Result Found:")
                # The "equation" is the statement of the binding affinity.
                # We output the key number here as requested.
                print(f"Binding Affinity (IC50) = {ic50_value} nM")
                
                # Step 5: Check which range the value falls into
                if 0.1 <= ic50_value <= 100:
                    print(f"This value falls in the range B: 0.1 - 100 nM")
                elif ic50_value < 0.1:
                    print(f"This value falls in the range A: < 0.1 nM")
                else:
                    print(f"This value is outside of the most common affinity ranges.")
                
                found_nM_value = True
                break # Stop after finding the first relevant result

        if not found_nM_value:
            print("No IC50 data in nM units was found in the database.")

    except Exception as e:
        print(f"An error occurred while connecting to the ChEMBL database: {e}")
        print("Please check your internet connection and ensure the client library is up to date.")

if __name__ == "__main__":
    find_binding_affinity()
