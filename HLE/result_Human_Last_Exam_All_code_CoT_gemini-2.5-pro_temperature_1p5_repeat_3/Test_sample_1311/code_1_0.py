import sys
import platform

# Check for Python version and library installation
try:
    if sys.version_info < (3, 6):
        print("This script requires Python 3.6 or newer.")
        sys.exit(1)
        
    from chembl_webresource_client.new_client import new_client
except ImportError:
    print("Required library 'chembl_webresource_client' not found.")
    print("Please install it by running the following command in your shell:")
    print("pip install chembl_webresource_client")
    sys.exit(1)
except Exception as e:
    print(f"An error occurred during import: {e}")
    sys.exit(1)

def find_binding_affinity():
    """
    Queries the ChEMBL database to find the binding affinity of a specific
    compound to CDK7 and determines which of the given ranges it falls into.
    """
    print("Querying ChEMBL database...")

    # Using pre-identified ChEMBL IDs for robustness and speed.
    # Compound: 7-dimethylphosphoryl-3-[2-[[(3~S)-6,6-dimethylpiperidin-3-yl]amino]-5-(trifluoromethyl)pyrimidin-4-yl]-1~H-indole-6-carbonitrile
    # This compound is also known as SY-1365.
    # Target: Cyclin-dependent kinase 7 (Homo sapiens)
    compound_chembl_id = 'CHEMBL3989895'
    target_chembl_id = 'CHEMBL277'
    
    try:
        # Initialize the activity client
        activity = new_client.activity

        # Filter for IC50 activity data for the specific compound-target pair
        res = activity.filter(
            molecule_chembl_id=compound_chembl_id,
            target_chembl_id=target_chembl_id,
            type='IC50',
            assay_type='B' # 'B' for Binding assays
        ).only(['standard_type', 'standard_relation', 'standard_value', 'standard_units'])

        if not res:
            print(f"No binding assay (IC50) data found for compound {compound_chembl_id} against target {target_chembl_id}.")
            return

        # Take the first and most relevant result
        activity_data = res[0]
        affinity_type = activity_data.get('standard_type')
        relation = activity_data.get('standard_relation')
        value_str = activity_data.get('standard_value')
        units = activity_data.get('standard_units')

        if not all([affinity_type, relation, value_str, units]):
            print("Retrieved data is incomplete. Cannot determine affinity.")
            return

        # Print the retrieved affinity value as an "equation"
        print("\nFound binding affinity from a published study:")
        print(f"Binding Affinity ({affinity_type}) {relation} {value_str} {units}")
        
        # Convert value to float for comparison
        value = float(value_str)

        # Determine which category the value falls into
        print("\nComparing this value to the provided answer choices:")
        answer_choice = "Unknown"
        if units == 'nM':
            if value < 0.1:
                answer_choice = "A"
            elif 0.1 <= value <= 100:
                answer_choice = "B"
        elif units == 'uM':
            if 0.1 <= value <= 100:
                answer_choice = "C"
        elif units == 'mM':
            if 0.1 <= value <= 100:
                answer_choice = "D"
            elif value > 100:
                answer_choice = "E"
        
        if answer_choice != "Unknown":
            print(f"The result {value} {units} falls into the range for choice {answer_choice}: '0.1 - 100 nM'.")
        else:
            print(f"Could not classify the result {value} {units} into the given choices.")

    except Exception as e:
        print(f"An error occurred while querying ChEMBL: {e}")
        print("Please check your internet connection and try again.")

if __name__ == "__main__":
    find_binding_affinity()