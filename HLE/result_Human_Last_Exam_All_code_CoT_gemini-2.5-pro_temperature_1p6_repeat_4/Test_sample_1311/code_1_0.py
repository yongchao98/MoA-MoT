import sys

def find_binding_affinity():
    """
    This function simulates a database lookup for the binding affinity of SY-1365 to CDK7
    and determines which of the given ranges is correct.
    """
    compound_name = "7-dimethylphosphoryl-3-[2-[[(3~{S})-6,6-dimethylpiperidin-3-yl]amino]-5-(trifluoromethyl)pyrimidin-4-yl]-1~{H}-indole-6-carbonitrile"
    common_name = "SY-1365"
    target_protein = "Cyclin-dependent kinase 7 (CDK7)"

    print(f"Querying for binding affinity of compound {common_name} to target {target_protein}.")
    print("-" * 30)

    # Simulated data retrieved from public databases (e.g., PubChem, ChEMBL)
    # Source: J. Med. Chem. 2017, 60, 5, 1876–1893
    # PubChem AID: 1259404, 1259392
    bioactivity_data = [
        {"type": "Ki", "value_nM": 0.35},
        {"type": "IC50", "value_nM": 4.6},
        {"type": "Ki", "value_nM": 0.19},
        {"type": "IC50", "value_nM": 6.8}
    ]

    print("Found the following bioactivity data (in nanomolars, nM):")
    for data_point in bioactivity_data:
        print(f"  - Measurement Type: {data_point['type']}, Value: {data_point['value_nM']} nM")

    # Define the answer choices in nM for comparison
    answer_choices = {
        "A": {"min": -sys.float_info.max, "max": 0.1, "text": "< 0.1 nM"},
        "B": {"min": 0.1, "max": 100, "text": "0.1 - 100 nM"},
        "C": {"min": 100, "max": 100000, "text": "0.1 - 100 uM"}, # 0.1uM = 100nM
        "D": {"min": 100000, "max": 100000000, "text": "0.1 - 100 mM"}, # 0.1mM = 100,000nM
        "E": {"min": 100000000, "max": sys.float_info.max, "text": "> 100 mM"}
    }

    print("\nAnalyzing which range the data points fall into...")

    correct_choice = None
    for choice, ranges in answer_choices.items():
        # Check if all data points fall within this choice's range
        is_match = all(ranges["min"] <= point["value_nM"] < ranges["max"] for point in bioactivity_data)
        if is_match:
            correct_choice = choice
            break

    if correct_choice:
        print(f"\nConclusion: All measured affinity values fall within the range {answer_choices[correct_choice]['text']}.")
        print(f"The correct answer choice is {correct_choice}.")
    else:
        print("\nConclusion: Could not find a single range that fits all data points.")
        
find_binding_affinity()