import sys
from chembl_webresource_client.new_client import new_client

def get_binding_affinity():
    """
    Connects to the ChEMBL database to find the binding affinity of a specific
    compound to a specific target protein.
    """
    # Define the unique identifiers for the compound and target protein.
    # Compound: 7-dimethylphosphoryl-3-[2-[[(3~{S})-6,6-dimethylpiperidin-3-yl]amino]-5-(trifluoromethyl)pyrimidin-4-yl]-1~{H}-indole-6-carbonitrile
    # (also known as SY-1365)
    compound_inchi_key = 'FFBYZBHUDGUPCR-CQSZACIVSA-N'
    
    # Target: Cyclin-dependent kinase 7 (CDK7)
    target_uniprot_id = 'P50613'

    print("Step 1: Connecting to ChEMBL database...")
    
    try:
        # Search for the molecule using its InChIKey to get its ChEMBL ID
        molecule = new_client.molecule.filter(molecule_structures__standard_inchi_key=compound_inchi_key).only(['molecule_chembl_id', 'pref_name'])
        if not molecule:
            print(f"Could not find molecule with InChIKey: {compound_inchi_key}")
            return
        compound_info = molecule[0]
        compound_chembl_id = compound_info['molecule_chembl_id']
        compound_name = compound_info.get('pref_name', 'N/A')

        # Search for the target using its UniProt ID to get its ChEMBL ID
        target = new_client.target.filter(target_components__accession=target_uniprot_id).only(['target_chembl_id', 'pref_name'])
        if not target:
            print(f"Could not find target with UniProt ID: {target_uniprot_id}")
            return
        # We are interested in the single protein target, not complexes
        target_info = [t for t in target if t['target_chembl_id'] == 'CHEMBL290'][0]
        target_chembl_id = target_info['target_chembl_id']
        target_name = target_info.get('pref_name', 'N/A')

        print(f"Step 2: Successfully identified compound '{compound_name}' ({compound_chembl_id})")
        print(f"Step 3: Successfully identified target '{target_name}' ({target_chembl_id})")
        print("-" * 40)
        
        # Search for activity data (IC50) for the compound-target pair
        print("Step 4: Querying for binding affinity data (IC50)...")
        activities = new_client.activity.filter(
            molecule_chembl_id=compound_chembl_id,
            target_chembl_id=target_chembl_id,
            standard_type="IC50",
            relation='='
        ).only(['standard_type', 'standard_value', 'standard_units'])

        if activities:
            print("\nFound the following IC50 values:")
            nM_values = []
            for act in activities:
                # Ensure the value is a number and units are present
                if act['standard_value'] and act['standard_units'] == 'nM':
                    value = float(act['standard_value'])
                    units = act['standard_units']
                    # Here we print the numbers from the "equation" (database record)
                    print(f"-> Equation (Measurement): IC50 = {value} {units}")
                    nM_values.append(value)
            
            if not nM_values:
                 print("\nNo IC50 values in nM units were found in the database.")
                 return

            # Determine the correct range based on the found values
            print("\nStep 5: Comparing to the answer choices...")
            print("A. < 0.1 nM")
            print("B. 0.1 - 100 nM")
            print("C. 0.1 - 100 uM")
            print("D. 0.1 - 100 mM")
            print("E. > 100 mM")

            # Check if all found values fall within a specific range
            if all(0.1 <= val <= 100 for val in nM_values):
                print("\nConclusion: All recorded values fall within the range of 0.1 - 100 nM.")
            else:
                print("\nConclusion: The recorded values fall into different ranges.")

        else:
            print("No IC50 binding affinity data was found for this specific compound-target pair.")

    except Exception as e:
        print(f"\nAn error occurred. Please ensure you have the 'chembl_webresource_client' library installed.", file=sys.stderr)
        print("You can install it with: pip install chembl_webresource_client", file=sys.stderr)
        print(f"Error details: {e}", file=sys.stderr)

if __name__ == "__main__":
    get_binding_affinity()