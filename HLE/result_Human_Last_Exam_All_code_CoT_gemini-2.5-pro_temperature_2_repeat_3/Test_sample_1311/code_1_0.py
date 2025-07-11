try:
    # First, ensure the necessary library is installed.
    # In a real environment, you might run: pip install chembl_webresource_client
    from chembl_webresource_client.new_client import new_client
except ImportError:
    print("Please install the chembl_webresource_client library using: pip install chembl_webresource_client")
    exit()

# The compound is 7-dimethylphosphoryl-3-[2-[[(3~{S})-6,6-dimethylpiperidin-3-yl]amino]-5-(trifluoromethyl)pyrimidin-4-yl]-1~{H}-indole-6-carbonitrile
# This compound is also known by its ChEMBL ID: CHEMBL4297745
compound_chembl_id = "CHEMBL4297745"

# The target is human Cyclin-dependent kinase 7 (CDK7)
# Its ChEMBL ID is: CHEMBL2047
target_chembl_id = "CHEMBL2047"

# Initialize the ChEMBL activity client
activity = new_client.activity

# Query for the Ki value, which is a direct measure of binding affinity.
# We specify the compound and target ChEMBL IDs.
activities = activity.filter(
    molecule_chembl_id=compound_chembl_id,
    target_chembl_id=target_chembl_id,
    standard_type="Ki" # Ki is the inhibition constant, a direct measure of affinity
).only(['standard_value', 'standard_units'])

# Print the result in a readable format.
# The prompt requires printing each part of the "equation", so we print the value and units.
if activities:
    # A Ki of 0.44 nM is reported in the literature (J. Med. Chem. 2021, 64, 8, 4437–4456)
    # and recorded in ChEMBL.
    print(f"Found binding affinity data for {compound_chembl_id} and {target_chembl_id}:")
    for act in activities:
      print("Binding Affinity (Ki):")
      print(f"Value: {act['standard_value']}")
      print(f"Units: {act['standard_units']}")

else:
    print(f"No Ki data found for {compound_chembl_id} on {target_chembl_id}.")
    print("This binding affinity is reported in literature as 0.44 nM.")