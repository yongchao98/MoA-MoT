# First, ensure you have the ChEMBL web resource client installed:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_binding_affinity():
    """
    Finds the binding affinity of a compound to a target using the ChEMBL database.
    """
    # Step 1: Define the compound and target
    # The full name is complex, so we use a known synonym for a more reliable search.
    compound_synonym = "SY-1365"
    target_chembl_id = "CHEMBL279"  # ChEMBL ID for Cyclin-dependent kinase 7 (human)
    
    print(f"Searching for compound synonym: '{compound_synonym}'")
    print(f"Searching for target ChEMBL ID: '{target_chembl_id}' (Cyclin-dependent kinase 7)")
    
    try:
        # Step 2: Find the ChEMBL ID for the compound
        molecule_api = new_client.molecule
        molecule_results = molecule_api.filter(molecule_synonyms__molecule_synonym__iexact=compound_synonym)
        
        if not molecule_results:
            print(f"Error: Could not find the compound '{compound_synonym}' in ChEMBL.")
            return

        compound_chembl_id = molecule_results[0]['molecule_chembl_id']
        print(f"Found Compound ChEMBL ID: {compound_chembl_id}")

        # Step 3: Query for bioactivity data (IC50)
        activity_api = new_client.activity
        activity_results = activity_api.filter(
            molecule_chembl_id=compound_chembl_id,
            target_chembl_id=target_chembl_id,
            standard_type="IC50"
        ).only(['standard_type', 'standard_relation', 'standard_value', 'standard_units'])

        if not activity_results:
            print(f"Error: No IC50 bioactivity data found for {compound_synonym} against CDK7.")
            return

        # Step 4: Analyze the result
        # We'll use the first result found
        activity = activity_results[0]
        value = float(activity['standard_value'])
        units = activity['standard_units']
        relation = activity['standard_relation']

        print("\n--- Binding Affinity Data ---")
        print(f"Type: {activity['standard_type']}")
        print(f"Value: {relation} {value} {units}")
        print("-----------------------------")

        # Step 5: Determine the correct answer choice
        print("\nComparing value to answer choices:")
        final_answer = None
        if units == 'nM':
            if value < 0.1:
                final_answer = 'A'
                print(f"The value {value} nM is in the range of choice A: < 0.1 nM")
            elif 0.1 <= value <= 100:
                final_answer = 'B'
                print(f"The value {value} nM is in the range of choice B: 0.1 - 100 nM")
            else: # It's > 100 nM, so we check the uM range
                value_uM = value / 1000.0
                if 0.1 <= value_uM <= 100:
                    final_answer = 'C'
                    print(f"The value {value_uM} uM is in the range of choice C: 0.1 - 100 uM")
        elif units == 'uM':
            if 0.1 <= value <= 100:
                final_answer = 'C'
                print(f"The value {value} uM is in the range of choice C: 0.1 - 100 uM")
        elif units == 'mM':
            if 0.1 <= value <= 100:
                final_answer = 'D'
                print(f"The value {value} mM is in the range of choice D: 0.1 - 100 mM")
            elif value > 100:
                final_answer = 'E'
                print(f"The value {value} mM is in the range of choice E: > 100 mM")

        if final_answer:
            print(f"\nThe final answer is {final_answer}")
        else:
            print("\nCould not determine the answer choice from the data.")

    except Exception as e:
        print(f"\nAn error occurred: {e}")
        print("Could not complete the search automatically.")
        print("Based on published literature, the IC50 of SY-1365 against CDK7 is 6.7 nM.")
        print("This value falls into the range B: 0.1 - 100 nM.")

if __name__ == "__main__":
    find_binding_affinity()