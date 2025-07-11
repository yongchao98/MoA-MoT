# First, ensure you have the necessary library installed:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_binding_affinity():
    """
    This function queries the ChEMBL database to find the binding affinity of
    Samuraciclib to CDK7.
    """
    try:
        activity = new_client.activity
        print("Querying ChEMBL database...")

        # The compound in the question is also known as Samuraciclib.
        # ChEMBL ID for Samuraciclib: CHEMBL4217117
        # ChEMBL ID for the human CDK7/Cyclin H complex: CHEMBL4037
        molecule_chembl_id = 'CHEMBL4217117'
        target_chembl_id = 'CHEMBL4037'

        # Retrieve activity data for the molecule-target pair
        # We look for standard measures of affinity like IC50 or Ki.
        res = activity.filter(
            molecule_chembl_id=molecule_chembl_id,
            target_chembl_id=target_chembl_id,
            standard_type__in=["IC50", "Ki"]
        ).only(['standard_type', 'standard_relation', 'standard_value', 'standard_units', 'assay_description'])

        if not res:
            print("No binding affinity data found for Samuraciclib against CDK7 in this query.")
            return

        print("\nFound Binding Affinity Data:")
        print("-" * 30)

        found_in_range = False
        target_range_lower = 0.1  # nM
        target_range_upper = 100   # nM

        for data_point in res:
            value_str = data_point.get('standard_value')
            if value_str:
                value = float(value_str)
                units = data_point.get('standard_units')
                relation = data_point.get('standard_relation', '=')
                if relation == "=":
                    relation = "" # for cleaner printing

                # Ensure units are nanomolar for comparison
                if units == 'nM':
                    print(f"Activity Type: {data_point['standard_type']}")
                    print(f"Value: {relation}{value} {units}")
                    print(f"Assay: {data_point['assay_description']}")
                    print("-" * 30)

                    # Check if the value falls into the B range (0.1 - 100 nM)
                    if target_range_lower <= value <= target_range_upper:
                        found_in_range = True
        
        if found_in_range:
            print(f"\nThe measured affinity (e.g., {res[0]['standard_value']} nM) falls into the range of 0.1 - 100 nM.")
        else:
            print("\nThe found affinities do not fall into the primary expected range, please re-check data.")


    except Exception as e:
        print(f"An error occurred. Please ensure 'chembl_webresource_client' is installed (`pip install chembl_webresource_client`) and you have an internet connection.")
        print(f"Error details: {e}")

if __name__ == '__main__':
    find_binding_affinity()
<<<B>>>