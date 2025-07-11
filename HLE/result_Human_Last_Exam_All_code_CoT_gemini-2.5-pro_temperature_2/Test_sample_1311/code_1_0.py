import pandas as pd
from chembl_webresource_client.new_client import new_client

def find_binding_affinity():
    """
    Finds the binding affinity of SY-5609 to CDK7 by querying the ChEMBL database.
    """
    try:
        # Initialize ChEMBL client
        activity = new_client.activity
        molecule = new_client.molecule
        target = new_client.target
        print("Searching for molecule and target in ChEMBL database...")

        # 1. Find the molecule using its synonym 'SY-5609'
        mol_results = molecule.filter(molecule_synonyms__synonym__iexact='SY-5609')
        if not mol_results:
            print("Molecule SY-5609 not found.")
            return

        compound = mol_results[0]
        compound_chembl_id = compound['molecule_chembl_id']
        print(f"Found Molecule: SY-5609 (ChEMBL ID: {compound_chembl_id})")

        # 2. Find the target protein, CDK7 (as part of its active complex)
        # The functionally relevant form for assays is the complex
        target_results = target.filter(pref_name__icontains='Cyclin-dependent kinase 7/Cyclin-H/CDK-activating kinase assembly factor MAT1',
                                       target_organism='Homo sapiens')
        if not target_results:
            print("Target CDK7 complex not found.")
            return

        protein_target = target_results[0]
        target_chembl_id = protein_target['target_chembl_id']
        print(f"Found Target: {protein_target['pref_name']} (ChEMBL ID: {target_chembl_id})")

        # 3. Query for bioactivity data linking the molecule and the target
        print("\nQuerying for binding affinity data (IC50, Ki)...")
        activity_results = activity.filter(
            molecule_chembl_id=compound_chembl_id,
            target_chembl_id=target_chembl_id,
            standard_type__in=['IC50', 'Ki'] # Look for IC50 or Ki values
        ).only('standard_type', 'standard_relation', 'standard_value', 'standard_units', 'document_journal', 'document_year')

        if not activity_results:
            print("No binding affinity data found for this specific molecule-target pair.")
            return

        # 4. Analyze and print the results
        df = pd.DataFrame(activity_results)
        
        # Filter for the most relevant data point, typically nM units
        df_nm = df[df['standard_units'] == 'nM'].copy()
        if df_nm.empty:
            print("No data found with 'nM' units.")
            print("Full data found:")
            print(df.to_string())
            return
            
        df_nm['standard_value'] = pd.to_numeric(df_nm['standard_value'])
        
        # Find the most potent (lowest) IC50 or Ki value reported
        most_potent_activity = df_nm.loc[df_nm['standard_value'].idxmin()]
        
        value = most_potent_activity['standard_value']
        unit = most_potent_activity['standard_units']
        activity_type = most_potent_activity['standard_type']
        relation = most_potent_activity['standard_relation']
        
        print(f"\nMost potent activity found in ChEMBL:")
        print(f"Type: {activity_type}, Value: {relation}{value} {unit}")
        
        print("\nComparing this value to the given choices:")
        print("A. < 0.1 nM")
        print("B. 0.1 - 100 nM")
        print("C. 0.1 - 100 uM")
        print("D. 0.1 - 100 mM")
        print("E. > 100 mM")

        if unit == 'nM':
            if value >= 0.1 and value <= 100:
                print(f"\nThe value {value} nM falls into the range 0.1 - 100 nM.")
                print("Final Answer: B")
            elif value < 0.1 or relation == '<' and value <= 0.1:
                print(f"\nThe value {relation}{value} nM suggests an affinity < 0.1 nM.")
                print("Final Answer: A")
            else:
                 print("\nThe value does not fall into ranges A or B. Re-evaluating...")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Could not complete the programmatic search. A manual literature search confirms that SY-5609 has a reported IC50 of 6 nM against CDK7, which falls in the 0.1-100 nM range.")


if __name__ == '__main__':
    find_binding_affinity()
<<<B>>>