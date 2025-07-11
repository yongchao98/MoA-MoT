import requests
import json

def get_binding_affinity():
    """
    Retrieves and prints the binding affinity of SY-1365 to human CDK7 from the ChEMBL database.
    """
    # ChEMBL IDs for the molecule (SY-1365) and the target (human CDK7)
    molecule_chembl_id = "CHEMBL3824364"
    target_chembl_id = "CHEMBL290"

    # Construct the API request URL to get activity data
    # We filter for standard_type 'IC50' and non-null standard_value
    url = (
        f"https://www.ebi.ac.uk/chembl/api/data/activity.json?"
        f"molecule_chembl_id={molecule_chembl_id}&"
        f"target_chembl_id={target_chembl_id}&"
        f"standard_type=IC50&"
        f"standard_value__isnull=false"
    )

    print(f"Querying ChEMBL for binding affinity data...")
    print(f"Compound: 7-dimethylphosphoryl-3-[2-[[(3S)-6,6-dimethylpiperidin-3-yl]amino]-5-(trifluoromethyl)pyrimidin-4-yl]-1H-indole-6-carbonitrile (SY-1365)")
    print(f"Target: Cyclin-dependent kinase 7 (CDK7)\n")

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)

        data = response.json()
        activities = data.get('activities', [])

        if not activities:
            print("No binding affinity data found for the specified compound and target.")
            return

        print("Found binding affinity data:")
        for activity in activities:
            activity_type = activity.get('standard_type', 'N/A')
            value = activity.get('standard_value', 'N/A')
            units = activity.get('standard_units', 'N/A')
            relation = activity.get('standard_relation', '')

            # In the final output, each number in the "equation" (measurement) is displayed.
            print(f"Type: {activity_type}, Value: {relation} {value} {units}")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while querying the ChEMBL API: {e}")
    except json.JSONDecodeError:
        print("Failed to parse the response from the API.")

if __name__ == "__main__":
    get_binding_affinity()