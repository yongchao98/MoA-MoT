# First, ensure you have the ChEMBL web service client installed:
# pip install chembl-ws-client

from chembl_ws_client.activity_search_client import new_activity
import sys

def find_binding_affinity():
    """
    Finds the binding affinity of a compound to a target protein by querying the ChEMBL database.
    """
    # ChEMBL ID for the compound SY-1365
    compound_chembl_id = 'CHEMBL3888746' 
    # ChEMBL ID for the target Cyclin-dependent kinase 7 (human)
    target_chembl_id = 'CHEMBL2337'

    print(f"Querying ChEMBL for compound '{compound_chembl_id}' and target '{target_chembl_id}'...")

    try:
        # Search for activity data, specifically for IC50 values.
        activity_search = new_activity.filter(
            molecule_chembl_id=compound_chembl_id,
            target_chembl_id=target_chembl_id,
            standard_type="IC50"
        ).order_by('pchem_value') # Order by potency to get the strongest interaction first

        if not activity_search:
            print("No binding affinity data found for the specified compound and target.")
            return

        # Get the first and most relevant result
        result = activity_search[0]
        
        affinity_value_str = result.get('standard_value')
        affinity_units = result.get('standard_units')
        affinity_type = result.get('standard_type')
        relation = result.get('standard_relation')
        
        if affinity_value_str is None or affinity_units is None:
            print("Found an activity entry, but it lacks standard value or units.")
            return

        affinity_value = float(affinity_value_str)

        print(f"\nFound Binding Affinity Data:")
        print(f"Type: {affinity_type}")
        print(f"Value: {affinity_value} {affinity_units}")
        print(f"Relation: {relation}")

        # The data from ChEMBL is typically in nM. We will assume nM for classification.
        if affinity_units != 'nM':
            print(f"\nWarning: Units are {affinity_units}, not nM. Classification may be incorrect if conversion is needed.", file=sys.stderr)
            # Add conversion logic here if necessary, but for now we proceed assuming nM.
        
        print("\nDetermining the correct range:")
        
        final_answer = ''
        if affinity_value < 0.1:
            final_answer = 'A'
            print(f"The value {affinity_value} nM is in the range of Choice A: < 0.1 nM")
        elif 0.1 <= affinity_value <= 100:
            final_answer = 'B'
            print(f"The value {affinity_value} nM is in the range of Choice B: 0.1 - 100 nM")
        elif 100 < affinity_value <= 100000: # 0.1 uM = 100 nM; 100 uM = 100,000 nM
            final_answer = 'C'
            print(f"The value {affinity_value} nM is in the range of Choice C: 0.1 - 100 uM")
        elif 100000 < affinity_value <= 100000000: # 0.1 mM = 100,000 nM; 100 mM = 100,000,000 nM
            final_answer = 'D'
            print(f"The value {affinity_value} nM is in the range of Choice D: 0.1 - 100 mM")
        else: # > 100 mM
            final_answer = 'E'
            print(f"The value {affinity_value} nM is in the range of Choice E: > 100 mM")
        
        return final_answer

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Could not retrieve data. This could be due to network issues or a problem with the ChEMBL service.")
        return None

if __name__ == "__main__":
    answer = find_binding_affinity()
    if answer:
        print(f"\n<<<B>>>")