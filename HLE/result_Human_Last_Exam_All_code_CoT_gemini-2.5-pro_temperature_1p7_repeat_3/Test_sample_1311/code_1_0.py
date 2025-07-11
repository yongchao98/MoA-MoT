try:
    # First, ensure you have the library installed:
    # pip install chembl_webresource_client
    
    from chembl_webresource_client.new_client import new_client

    # Step 1: Initialize the client for accessing the 'activity' data
    activity_client = new_client.activity

    # Step 2: Define the compound and target identifiers.
    # The compound is SY-1365, which corresponds to ChEMBL ID CHEMBL3989823.
    # The target is human CDK7, which corresponds to ChEMBL ID CHEMBL2147.
    COMPOUND_ID = 'CHEMBL3989823'
    TARGET_ID = 'CHEMBL2147'

    print(f"Querying ChEMBL for the binding affinity of compound {COMPOUND_ID}")
    print(f"against target {TARGET_ID} (Human CDK7)...")

    # Step 3: Query for IC50 values linking the compound and target
    activities = activity_client.filter(
        molecule_chembl_id=COMPOUND_ID,
        target_chembl_id=TARGET_ID,
        standard_type="IC50"  # IC50 is a standard measure of a drug's potency
    )

    # Step 4: Find the most relevant activity data and print it.
    # We will look for the specific assay against the CDK7/cyclin H/MAT1 complex.
    found_activity = False
    if activities:
        for activity in activities:
            # Check for the specific, well-characterized assay result
            if activity['assay_description'] and "CDK7/cyclin H/MAT1" in activity['assay_description']:
                value = activity['standard_value']
                units = activity['standard_units']
                
                print("\n--- Found Result ---")
                print(f"Binding Affinity (IC50) = {value} {units}")
                
                # Check which answer choice this falls into
                if units == 'nM' and 0.1 <= float(value) <= 100:
                    print("This value falls within the range B: 0.1 - 100 nM")
                
                found_activity = True
                break # Stop after finding the first relevant result
    
    if not found_activity:
        print("\nCould not find a specific IC50 value for this compound-target pair in the database.")

except ImportError:
    print("The 'chembl_webresource_client' library is not installed.")
    print("Please install it by running: pip install chembl_webresource_client")
except Exception as e:
    print(f"An error occurred: {e}")
    print("This may be due to a network issue or a change in the ChEMBL API.")
