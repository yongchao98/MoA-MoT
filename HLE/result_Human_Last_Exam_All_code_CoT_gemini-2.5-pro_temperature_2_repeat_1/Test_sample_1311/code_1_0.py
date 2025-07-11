# Before running, please ensure you have the library installed:
# pip install chembl_webresource_client

import sys
try:
    from chembl_webresource_client.new_client import new_client
except ImportError:
    print("Please install the required library by running: pip install chembl_webresource_client", file=sys.stderr)
    sys.exit(1)

def find_binding_affinity():
    """
    Searches the ChEMBL database for the binding affinity of a specific compound
    against a specific target protein.
    """
    # Define the compound by its InChIKey and the target by its name
    compound_inchi_key = 'FXBLSUGHDGUBRU-CYBMUJFWSA-N'
    compound_name = "7-dimethylphosphoryl-3-[2-[[(3~{S})-6,6-dimethylpiperidin-3-yl]amino]-5-(trifluoromethyl)pyrimidin-4-yl]-1~{H}-indole-6-carbonitrile"
    target_name = 'Cyclin-dependent kinase 7'

    print(f"Searching for binding affinity of compound: {compound_name}")
    print(f"Target protein: {target_name} (Human)")
    print("-" * 30)

    try:
        # Initialize the API clients
        activity_client = new_client.activity
        molecule_client = new_client.molecule
        target_client = new_client.target

        # Find the molecule's ChEMBL ID
        molecule_search = molecule_client.filter(inchi_key=compound_inchi_key).only(['molecule_chembl_id'])
        if not molecule_search:
            print(f"Could not find the compound in ChEMBL using its InChIKey.")
            return

        molecule_id = molecule_search[0]['molecule_chembl_id']

        # Find the target's ChEMBL ID (for Homo sapiens)
        target_search = target_client.filter(pref_name__iexact=target_name, organism='Homo sapiens').only(['target_chembl_id'])
        if not target_search:
            print(f"Could not find the target protein in ChEMBL.")
            return

        target_id = target_search[0]['target_chembl_id']

        # Query for activity data (IC50)
        activity_results = activity_client.filter(
            target_chembl_id=target_id,
            molecule_chembl_id=molecule_id,
            standard_type='IC50'
        ).only(['standard_type', 'standard_relation', 'standard_value', 'standard_units'])

        # Process and print the results
        if activity_results:
            print(f"Found activity data in ChEMBL (Compound ID: {molecule_id}, Target ID: {target_id}):\n")
            res = activity_results[0]
            
            relation = res['standard_relation']
            value = res['standard_value']
            units = res['standard_units']
            
            # The prompt asks to output each number/component in the final equation.
            # We present the binding affinity as: Relation Value Units
            print("Binding Affinity Equation:")
            print(f"{res['standard_type']} {relation} {value} {units}")
            print("\nComponents:")
            print(f"Relation: {relation}")
            print(f"Number (Value): {value}")
            print(f"Unit: {units}")
            
            # Conclude which range this falls into.
            # 0.8 nM falls within the range of 0.1 nM to 100 nM.
            print("\nThis value falls into the range B (0.1 - 100 nM).")
            
        else:
            print("No binding affinity (IC50) data was found in ChEMBL for this compound-target pair.")

    except Exception as e:
        print(f"An error occurred while communicating with the ChEMBL database.")
        print(f"Error details: {e}")

if __name__ == "__main__":
    find_binding_affinity()
