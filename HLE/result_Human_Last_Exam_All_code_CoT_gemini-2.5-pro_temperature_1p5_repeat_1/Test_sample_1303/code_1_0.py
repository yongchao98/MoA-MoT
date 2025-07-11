# First, please make sure you have the ChEMBL web resource client installed.
# If not, run the following command in your terminal:
# pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client

def find_interacting_molecules(target_id):
    """
    Finds and prints small molecules that interact with a given ChEMBL target ID.

    This function queries the ChEMBL database for activities associated with the
    target ID, filtering for potent interactions (pChEMBL >= 5), and then prints
    the details of each interacting molecule.
    """
    try:
        activity = new_client.activity
        # Query for activities with a pChEMBL value of 5 or greater.
        # This corresponds to a potency of 10 uM or less (e.g., IC50, Ki).
        # We also filter for common standard types to get quantitative data.
        results = activity.filter(
            target_chembl_id=target_id, 
            pchembl_value__gte=5,
            standard_type__in=['IC50', 'Ki', 'Kd', 'EC50']
        ).only([
            'molecule_chembl_id', 
            'standard_type', 
            'standard_value', 
            'standard_units'
        ])

        if not results:
            print(f"No interacting small molecules found for ChEMBL ID: {target_id} with the specified criteria (pChEMBL >= 5).")
            return

        print(f"Found interacting small molecules for target {target_id}:")
        
        # Use a set to avoid printing duplicate molecule IDs if they have multiple activity entries
        processed_molecules = set()
        
        for interaction in results:
            molecule_id = interaction.get('molecule_chembl_id')
            if molecule_id and molecule_id not in processed_molecules:
                activity_type = interaction.get('standard_type', 'N/A')
                value = interaction.get('standard_value', 'N/A')
                units = interaction.get('standard_units', '')
                
                print(f"- Molecule: {molecule_id}, Interaction Type: {activity_type}, Value: {value} {units}")
                processed_molecules.add(molecule_id)

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please check your internet connection and ensure the 'chembl_webresource_client' library is correctly installed.")

if __name__ == "__main__":
    # The ChEMBL ID for the target of interest
    target_chembl_id = "CHEMBL4803817"
    find_interacting_molecules(target_chembl_id)