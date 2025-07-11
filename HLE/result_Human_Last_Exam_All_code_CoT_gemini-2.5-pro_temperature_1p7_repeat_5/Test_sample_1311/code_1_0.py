import sys
import warnings

# Suppress all warnings from the ChEMBL client for a cleaner output
warnings.filterwarnings("ignore")

# Step 1: Install and import the necessary ChEMBL client library
try:
    from chembl_webresource_client.new_client import new_client
except ImportError:
    print("The 'chembl_webresource_client' library is not installed.")
    print("Please run the following command to install it:")
    print("pip install chembl_webresource_client")
    sys.exit(1)

# Step 2: Define identifiers for the molecule and target protein
# InChIKey is a standard, unique identifier for a chemical structure.
# This key corresponds to '7-dimethylphosphoryl-3-[2-[[(3~{S})-6,6-dimethylpiperidin-3-yl]amino]-5-(trifluoromethyl)pyrimidin-4-yl]-1~{H}-indole-6-carbonitrile'
compound_inchi_key = "QGYSGGPWFJPCGV-QTGZJDMSSA-N"
target_name = "Cyclin-dependent kinase 7"

# Initialize ChEMBL API clients
molecule_client = new_client.molecule
target_client = new_client.target
activity_client = new_client.activity

# Step 3: Find the molecule and target in ChEMBL
try:
    # Find the molecule using its InChIKey for accuracy
    molecule_results = molecule_client.filter(molecule_structures__standard_inchi_key=compound_inchi_key)
    if not molecule_results:
        print(f"Error: Could not find the specified molecule in the ChEMBL database.")
        sys.exit(1)
    molecule_data = molecule_results[0]
    molecule_chembl_id = molecule_data['molecule_chembl_id']
    molecule_name = molecule_data.get('pref_name', "the specified compound")

    # Find the human version of the target protein
    target_results = target_client.search(target_name).filter(organism="Homo sapiens")
    if not target_results:
        print(f"Error: Could not find the human target '{target_name}' in the ChEMBL database.")
        sys.exit(1)
    target_data = target_results[0]
    target_chembl_id = target_data['target_chembl_id']
    target_pref_name = target_data['pref_name']

    # Step 4: Query for binding affinity data linking the molecule and target
    print(f"Searching for binding affinity between '{molecule_name}' and '{target_pref_name}'...")
    
    # Filter for common affinity/potency measurements (IC50, Ki, Kd)
    activity_results = activity_client.filter(
        molecule_chembl_id=molecule_chembl_id,
        target_chembl_id=target_chembl_id,
        standard_type__in=["IC50", "Ki", "Kd"]
    ).only(['standard_type', 'standard_relation', 'standard_value', 'standard_units'])

    if not activity_results:
        print("No binding affinity data found in ChEMBL for this specific compound-target pair.")
        sys.exit(1)

    # Step 5: Analyze the data and present the conclusion
    # We'll use the first valid result found
    activity_record = activity_results[0]
    act_type = activity_record['standard_type']
    relation = activity_record['standard_relation']
    value_str = activity_record['standard_value']
    units = activity_record['standard_units']

    if not all([act_type, relation, value_str, units]):
         print("Found an activity record, but it is missing complete data for analysis.")
         sys.exit(1)

    value = float(value_str)

    print("\n--- Bioactivity Data Found ---")
    print(f"The binding affinity is expressed as:")
    # As requested, output each number in the final equation
    print(f"{act_type} {relation} {value} {units}")
    print("----------------------------\n")

    # Step 6: Compare with answer choices to find the correct range
    final_answer = ""
    final_answer_text = ""

    if units == 'nM': # Nanomolar
        if value < 0.1:
            final_answer = 'A'
            final_answer_text = "< 0.1 nM"
        elif 0.1 <= value <= 100:
            final_answer = 'B'
            final_answer_text = "0.1 - 100 nM"
        else: # > 100 nM, check if it falls into uM range
            value_in_uM = value / 1000
            if 0.1 <= value_in_uM <= 100:
                final_answer = 'C'
                final_answer_text = "0.1 - 100 uM"
    elif units == 'uM': # Micromolar
        if 0.1 <= value <= 100:
            final_answer = 'C'
            final_answer_text = "0.1 - 100 uM"
        else: # Check if it falls in other ranges
            value_in_nM = value * 1000
            if value_in_nM < 0.1:
                final_answer = 'A'
                final_answer_text = "< 0.1 nM"
            elif 0.1 <= value_in_nM <= 100:
                final_answer = 'B'
                final_answer_text = "0.1 - 100 nM"
    
    if final_answer:
        print(f"Conclusion: The value {value} {units} falls into choice {final_answer} ({final_answer_text}).")
    else:
        print(f"The value {value} {units} does not fit into the provided answer choices (A, B, C).")

except Exception as e:
    print(f"An error occurred during the process: {e}")
    print("This may be due to a temporary issue with the ChEMBL web service or a network problem.")
