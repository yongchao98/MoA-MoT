import requests

def find_binding_affinity():
    """
    Finds the binding affinity of Samuraciclib to CDK7 by querying the ChEMBL database.
    """
    # ChEMBL IDs for the molecule and target
    molecule_chembl_id = "CHEMBL3988950"
    molecule_name = "Samuraciclib"
    target_chembl_id = "CHEMBL301"
    target_name = "Cyclin-dependent kinase 7"

    # API endpoint to fetch IC50 activity data
    url = (
        f"https://www.ebi.ac.uk/chembl/api/data/activity.json?"
        f"molecule_chembl_id={molecule_chembl_id}&"
        f"target_chembl_id={target_chembl_id}&"
        f"standard_type=IC50"
    )

    print("Querying ChEMBL database for binding affinity...")
    print(f"Compound: {molecule_name} ({molecule_chembl_id})")
    print(f"Target: {target_name} ({target_chembl_id})")
    print("-" * 20)

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()  # Will raise an HTTPError for bad responses
        data = response.json()

        if not data or not data.get('activities'):
            print("No bioactivity data found for the specified molecule and target.")
            return

        # Use the first activity entry as a representative value
        activity = data['activities'][0]
        value_str = activity.get('standard_value')
        units = activity.get('standard_units')
        relation = activity.get('standard_relation', '')

        if not value_str or not units:
            print("Found activity data, but it's missing standard value or units.")
            return

        value = float(value_str)

        # Print the equation/result in a clear format
        print(f"Found Activity: IC50 {relation} {value} {units}")
        print("-" * 20)
        
        # Determine which category the value falls into
        print("Comparing this value to the answer choices:")
        choices = {
            "A": "< 0.1 nM",
            "B": "0.1 - 100 nM",
            "C": "0.1 - 100 uM",
            "D": "0.1 - 100 mM",
            "E": "> 100 mM"
        }
        for key, text in choices.items():
            print(f"{key}. {text}")

        answer = None
        if units == "nM":
            if value < 0.1:
                answer = "A"
            elif 0.1 <= value <= 100:
                answer = "B"
            elif 100 < value < 100000: # 0.1-100 uM
                answer = "C"
            
        elif units == "uM":
            if 0.1 <= value <= 100:
                answer = "C"

        if answer:
            print(f"\nThe value {value} {units} falls into the range '{choices[answer]}'.")
            print(f"The correct answer is {answer}.")
        else:
            print(f"\nCould not classify the value {value} {units} into the given choices.")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while contacting the ChEMBL API: {e}")
    except (KeyError, IndexError, ValueError) as e:
        print(f"An error occurred while parsing the data: {e}")

if __name__ == '__main__':
    find_binding_affinity()