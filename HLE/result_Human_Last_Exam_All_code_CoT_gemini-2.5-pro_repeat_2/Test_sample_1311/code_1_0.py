# First, ensure you have the ChEMBL web resource client installed.
# If not, run: pip install chembl_webresource_client
try:
    from chembl_webresource_client.new_client import new_client
except ImportError:
    print("Please install the ChEMBL client first by running: pip install chembl_webresource_client")
    exit()

# The compound is Samuraciclib (CHEMBL4297690)
# The target is Cyclin-dependent kinase 7 (human, CHEMBL301)
COMPOUND_ID = 'CHEMBL4297690'
TARGET_ID = 'CHEMBL301'
ACTIVITY_TYPE = 'IC50'
ACTIVITY_UNIT = 'nM'

print(f"Querying ChEMBL database for binding affinity...")
print(f"Compound (Samuraciclib): {COMPOUND_ID}")
print(f"Target (CDK7): {TARGET_ID}")
print("-" * 30)

try:
    # Set up the activity client
    activity = new_client.activity

    # Filter for the specific activity data
    res = activity.filter(
        molecule_chembl_id=COMPOUND_ID,
        target_chembl_id=TARGET_ID,
        standard_type=ACTIVITY_TYPE,
        standard_units=ACTIVITY_UNIT
    ).only(['standard_value', 'standard_relation'])

    if not res:
        print("No definitive IC50 data in nM found for this compound-target pair.")
    else:
        # We will take the first result as a representative value
        # ChEMBL can contain multiple results from different assays
        result = res[0]
        ic50_value = result['standard_value']
        relation = result['standard_relation']
        
        print(f"Found binding affinity data:")
        print(f"The reported IC50 value is: {relation} {ic50_value} {ACTIVITY_UNIT}")
        
        # The prompt asks to output each number in the final equation.
        # Here, the final result is the value itself.
        print("\nFinal Result Number:")
        print(ic50_value)

except Exception as e:
    print(f"An error occurred: {e}")
    print("Could not retrieve data from ChEMBL. Please check your internet connection and if the service is available.")
