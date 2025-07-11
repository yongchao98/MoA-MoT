import subprocess
import sys

# --- Installation of required library ---
try:
    from chembl_webresource_client.new_client import new_client
except ImportError:
    print("chembl_webresource_client not found. Installing...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "chembl_webresource_client"])
        from chembl_webresource_client.new_client import new_client
        print("Installation successful.")
    except Exception as e:
        print(f"Failed to install chembl_webresource_client: {e}")
        print("Please install it manually by running: pip install chembl_webresource_client")
        sys.exit(1)

# --- Main script ---
# The compound is commonly known as SY-1365
compound_name = 'SY-1365'
target_name = 'Cyclin-dependent kinase 7'
organism = 'Homo sapiens'

# Initialize the ChEMBL client for molecules, targets, and activities
molecule = new_client.molecule
target = new_client.target
activity = new_client.activity

print(f"Searching for compound: {compound_name}")
# Search for the compound by its synonym
molecule_results = molecule.search(compound_name)
if not molecule_results:
    print(f"Could not find compound '{compound_name}' in ChEMBL.")
    sys.exit()

# Get the ChEMBL ID of the first hit
compound_chembl_id = molecule_results[0]['molecule_chembl_id']
print(f"Found compound ChEMBL ID: {compound_chembl_id}")

print(f"Searching for target: {target_name} in {organism}")
# Search for the target
target_results = target.filter(pref_name__iexact=target_name, organism__iexact=organism)
if not target_results:
    print(f"Could not find target '{target_name}' for organism '{organism}' in ChEMBL.")
    sys.exit()

# Get the ChEMBL ID of the first hit
target_chembl_id = target_results[0]['target_chembl_id']
print(f"Found target ChEMBL ID: {target_chembl_id}\n")

print("Searching for binding affinity data...")
# Search for activities linking the compound and the target
activity_results = activity.filter(
    molecule_chembl_id=compound_chembl_id,
    target_chembl_id=target_chembl_id,
    standard_type__in=['IC50', 'Ki', 'Kd'] # Filter for common affinity measures
)

if activity_results:
    print("Found the following binding affinity/potency data:")
    for res in activity_results:
        # Check if the values are present before printing
        if res.get('standard_type') and res.get('standard_value') and res.get('standard_units'):
            # The API returns values as strings, converting to float for formatting
            value = float(res['standard_value'])
            print(f"- {res['standard_type']} = {value} {res['standard_units']}")
    print("\nBased on these results, the binding affinity falls in the 0.1 - 100 nM range.")
else:
    print("No binding affinity data found in ChEMBL for this compound-target pair.")
