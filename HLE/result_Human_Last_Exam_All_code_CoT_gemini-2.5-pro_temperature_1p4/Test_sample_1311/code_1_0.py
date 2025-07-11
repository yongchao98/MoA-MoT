import sys

try:
    # We use the ChEMBL web resource client to query the database.
    from chembl_webresource_client.new_client import new_client
except ImportError:
    print("Error: The 'chembl_webresource_client' library is not installed.")
    print("Please install it by running: pip install chembl_webresource_client")
    sys.exit(1)

# Step 1: Define identifiers for the compound and the target protein.
# The long chemical name corresponds to the investigational drug SY-5609.
compound_synonym = 'SY-5609'
# The UniProt accession ID for human Cyclin-dependent kinase 7 (CDK7) is P50613.
target_uniprot_id = 'P50613'

print(f"Searching for compound: {compound_synonym}")
print(f"Searching for target: Human CDK7 (UniProt ID: {target_uniprot_id})\n")

# Step 2: Initialize clients to access different parts of the ChEMBL database.
activity_client = new_client.activity
molecule_client = new_client.molecule
target_client = new_client.target

# Step 3: Find the ChEMBL ID for the molecule using its synonym.
mol_results = molecule_client.filter(molecule_synonyms__molecule_synonym__iexact=compound_synonym)
if not mol_results:
    print(f"Could not find a ChEMBL entry for molecule '{compound_synonym}'")
    sys.exit(1)
molecule_id = mol_results[0]['molecule_chembl_id']
print(f"Found Molecule ID: {molecule_id}")

# Step 4: Find the ChEMBL ID for the target using its UniProt accession ID.
target_results = target_client.filter(target_components__accession=target_uniprot_id, organism='Homo sapiens')
if not target_results:
    print(f"Could not find a ChEMBL entry for target with UniProt ID '{target_uniprot_id}'")
    sys.exit(1)
target_id = target_results[0]['target_chembl_id']
print(f"Found Target ID: {target_id}\n")

# Step 5: Query for activity data linking the molecule and target.
# We look for standard binding affinity measures like IC50 or Ki.
print("Querying for binding affinity data...")
activity_results = activity_client.filter(
    molecule_chembl_id=molecule_id,
    target_chembl_id=target_id,
    standard_type__in=["IC50", "Ki"]
).only(['standard_type', 'standard_relation', 'standard_value', 'standard_units'])

if not activity_results:
    print("No binding affinity data (IC50, Ki) was found for this compound-target pair.")
else:
    # Step 6: Print the most relevant results.
    # We prioritize exact measurements ('=') in nanomolar units ('nM').
    found_exact_match = False
    for record in activity_results:
        if record['standard_units'] == 'nM' and record['standard_relation'] == '=' and record['standard_value'] is not None:
            affinity_type = record['standard_type']
            affinity_value = record['standard_value']
            affinity_units = record['standard_units']
            
            # Print each component of the final result as requested
            print("\n--- Binding Affinity Found ---")
            print(f"Type: {affinity_type}")
            print(f"Value: {affinity_value}")
            print(f"Units: {affinity_units}")
            print(f"Final Equation: {affinity_type} = {affinity_value} {affinity_units}")
            found_exact_match = True
            # We can break after finding a good result, as they are often duplicated.
            break
            
    if not found_exact_match:
        print("\nCould not find an exact (relation '=') binding affinity in nM.")
        print("Showing all available records:")
        for record in activity_results:
            print(record)
