# First, ensure you have the ChEMBL client library installed:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_binding_affinity():
    """
    Finds the binding affinity of Zotiraciclib to CDK7 using the ChEMBL database.
    """
    try:
        # Define ChEMBL IDs for the compound (Zotiraciclib) and target (CDK7)
        compound_chembl_id = 'CHEMBL3835086'
        target_chembl_id = 'CHEMBL301'

        # Initialize the ChEMBL activity client
        activity = new_client.activity

        # Query for binding data (Ki or IC50) for the specified compound and target
        res = activity.filter(
            molecule_chembl_id=compound_chembl_id,
            target_chembl_id=target_chembl_id
        ).filter(
            standard_type__in=["Ki", "IC50"]  # Prioritize Ki, but accept IC50
        ).only(
            ['standard_type', 'standard_value', 'standard_units', 'standard_relation']
        )

        if not res:
            print("No binding affinity data found in ChEMBL for this compound-target pair.")
            return

        # Find the most relevant data point (preferring Ki)
        best_data = None
        for r in res:
            if r['standard_type'] == 'Ki' and r['standard_value'] is not None:
                best_data = r
                break # Found Ki, so stop looking
        
        # If no Ki found, take the first available record
        if not best_data and res:
            best_data = res[0]
            
        if not best_data or best_data['standard_value'] is None:
            print("Found assay entries, but they lack specific standard values.")
            return

        affinity_type = best_data['standard_type']
        affinity_value_str = best_data['standard_value']
        affinity_units = best_data['standard_units']
        affinity_relation = best_data['standard_relation']
        
        # Ensure units are 'nM' for comparison
        if affinity_units != 'nM':
            print(f"Data found in units of {affinity_units}, which cannot be directly compared.")
            return

        affinity_value = float(affinity_value_str)
        
        print(f"Found binding affinity from ChEMBL database:")
        print(f"- Type: {affinity_type}")
        print(f"- Value: {affinity_relation} {affinity_value} {affinity_units}")
        print("-" * 20)

        # Determine which category the value falls into
        if affinity_value < 0.1:
            print(f"The value {affinity_value} nM falls into choice A (< 0.1 nM).")
            print(f"Final Equation: {affinity_value} < 0.1")
        elif 0.1 <= affinity_value <= 100:
            print(f"The value {affinity_value} nM falls into choice B (0.1 - 100 nM).")
            print(f"Final Equation: 0.1 <= {affinity_value} <= 100")
        elif 100 < affinity_value <= 100000: # 0.1 uM = 100 nM, 100 uM = 100,000 nM
            print(f"The value {affinity_value} nM falls into choice C (0.1 - 100 uM).")
            # Equation for this would be based on uM, not clear how to represent mixed units.
        else:
            print("The value does not fit within the common nM or uM ranges specified.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have the 'chembl_webresource_client' library installed (`pip install chembl_webresource_client`) and an active internet connection.")

# Run the function
find_binding_affinity()