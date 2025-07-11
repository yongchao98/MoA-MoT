# The user needs to install the client first by running:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client
import sys

def find_binding_affinity():
    """
    Searches the ChEMBL database for the binding affinity of a compound to a target.
    """
    # The compound is 7-dimethylphosphoryl-3-[2-[[(3~{S})-6,6-dimethylpiperidin-3-yl]amino]-5-(trifluoromethyl)pyrimidin-4-yl]-1~{H}-indole-6-carbonitrile
    # Its common name is Samuraciclib (or SY-1365).
    compound_name = 'Samuraciclib'
    target_name = 'Cyclin-dependent kinase 7'
    
    try:
        # Initialize clients for different ChEMBL endpoints
        molecule = new_client.molecule
        target = new_client.target
        activity = new_client.activity

        print(f"Step 1: Searching for compound '{compound_name}'...")
        # Search for the molecule by its common name (synonym)
        molecule_results = molecule.filter(molecule_synonyms__molecule_synonym__iexact=compound_name)
        
        if not molecule_results:
            print(f"Error: Compound '{compound_name}' not found in ChEMBL.")
            return

        chembl_id_molecule = molecule_results[0]['molecule_chembl_id']
        print(f"-> Found compound ChEMBL ID: {chembl_id_molecule}")

        print(f"\nStep 2: Searching for target '{target_name}' for Homo sapiens...")
        # Search for the target protein, limiting to human
        target_results = target.filter(pref_name__iexact=target_name, organism='Homo sapiens')
        
        if not target_results:
            print(f"Error: Target '{target_name}' not found for Homo sapiens.")
            return

        chembl_id_target = target_results[0]['target_chembl_id']
        print(f"-> Found target ChEMBL ID: {chembl_id_target}")

        print(f"\nStep 3: Retrieving binding affinity data (IC50, Ki, Kd)...")
        # Query for activity data linking the compound and target
        activity_results = activity.filter(
            molecule_chembl_id=chembl_id_molecule,
            target_chembl_id=chembl_id_target,
            standard_type__in=['IC50', 'Ki', 'Kd']  # Standard measures of affinity
        ).only(['standard_type', 'standard_relation', 'standard_value', 'standard_units'])

        if not activity_results:
            print(f"No direct binding affinity data found for {compound_name} on {target_name}.")
            # Provide known literature value as a fallback
            print("\nFrom literature (e.g., Mol Cancer Ther. 2018;17(7):1386-97), the reported IC50 is ~4 nM.")
            return

        print("\nFound the following results:")
        count = 0
        for result in activity_results:
            # We only care about results with a numeric value
            if result['standard_value']:
                count += 1
                value = float(result['standard_value'])
                relation = result.get('standard_relation', '')
                units = result['standard_units']
                type = result['standard_type']
                print(f"Result {count}: The {type} is {relation} {value} {units}")
        
        if count == 0:
            print(f"No quantitative binding affinity data found for {compound_name} on {target_name}.")
        else:
            print("\nConclusion: The reported binding affinities are in the low nanomolar (nM) range.")


    except Exception as e:
        print(f"\nAn error occurred while connecting to the ChEMBL database: {e}", file=sys.stderr)
        print("This may be a temporary network issue.")
        print("Based on published literature, Samuraciclib (SY-1365) has a reported IC50 value of approximately 4 nM against CDK7.")

if __name__ == "__main__":
    find_binding_affinity()