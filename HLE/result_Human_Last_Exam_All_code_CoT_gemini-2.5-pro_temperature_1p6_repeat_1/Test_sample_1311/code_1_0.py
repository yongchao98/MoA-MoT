# First, ensure you have the ChEMBL client library installed:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_binding_affinity():
    """
    Finds the binding affinity of a compound to a target protein
    using the ChEMBL database.
    """
    compound_synonym = 'SY-5609'
    target_name = 'Cyclin-dependent kinase 7'
    target_organism = 'Homo sapiens'
    
    # Set up clients for ChEMBL's molecule, target, and activity resources
    molecule = new_client.molecule
    target = new_client.target
    activity = new_client.activity

    try:
        # Step 1: Find the compound's ChEMBL ID using its synonym
        print(f"Searching for compound with synonym: '{compound_synonym}'...")
        mol_results = molecule.filter(synonyms__icontains=compound_synonym)
        if not mol_results:
            print(f"Error: Could not find compound '{compound_synonym}' in ChEMBL.")
            return
        
        compound_chembl_id = mol_results[0]['molecule_chembl_id']
        print(f"Found Compound ID: {compound_chembl_id}")

        # Step 2: Find the target's ChEMBL ID
        print(f"Searching for target: '{target_name}' in organism '{target_organism}'...")
        target_results = target.filter(pref_name__iexact=target_name, organism=target_organism)
        if not target_results:
            print(f"Error: Could not find target '{target_name}' for '{target_organism}'.")
            return
            
        target_chembl_id = target_results[0]['target_chembl_id']
        print(f"Found Target ID: {target_chembl_id}")

        # Step 3: Query for activity between the compound and target
        print(f"Querying for activity between {compound_chembl_id} and {target_chembl_id}...")
        activity_results = activity.filter(
            molecule_chembl_id=compound_chembl_id,
            target_chembl_id=target_chembl_id,
            standard_type__in=['Ki', 'IC50'] # Look for Ki or IC50 values
        ).order_by('standard_value') # Get the most potent value first
        
        if not activity_results:
            print("Error: No binding affinity data (Ki or IC50) found for this pair.")
            return

        # Step 4: Analyze and report the result
        # We will use the first result, which is the most potent one found
        best_result = activity_results[0]
        
        affinity_type = best_result['standard_type']
        relation = best_result['standard_relation']
        value = best_result['standard_value']
        units = best_result['standard_units']
        
        if value is None or units is None:
            print("Found an activity entry, but the value or units are missing.")
            return
        
        print("\n--- Found Binding Affinity Data ---")
        # The prompt requested to output each number in the final equation.
        # Here we print the components of the measurement: Type, Relation, Value, and Unit.
        print(f"Measurement Type: {affinity_type}")
        print(f"Relation: {relation}")
        print(f"Value: {value}")
        print(f"Units: {units}")
        print("-----------------------------------")
        
        print(f"\nThe reported binding affinity is: {affinity_type} {relation} {value} {units}.")

        # Determine which range the value falls into
        val_float = float(value)
        
        # We are looking for values in the range of 0.1 to 100 nM
        if units == 'nM' and 0.1 <= val_float <= 100:
            final_answer = "B. 0.1 - 100 nM"
        elif units == 'nM' and val_float < 0.1:
            final_answer = "A. < 0.1 nM"
        # Convert uM to nM for comparison (1 uM = 1000 nM)
        elif units == 'uM' and 0.1 <= val_float <= 100:
            final_answer = "C. 0.1 - 100 uM"
        else:
            final_answer = "Value is outside the primary ranges, please check the data."
            
        print(f"\nThis value falls into the category: {final_answer}")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("This could be due to a network issue or a problem with the ChEMBL service.")

if __name__ == "__main__":
    find_binding_affinity()

<<<B>>>