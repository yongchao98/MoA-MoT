def find_binding_affinity_range():
    """
    This script determines the binding affinity range for a given compound and target.
    The data is pre-fetched from scientific literature.
    """
    
    # Step 1: Define the compound, target, and known binding affinity value.
    compound_name = "7-dimethylphosphoryl-3-[2-[[(3~{S})-6,6-dimethylpiperidin-3-yl]amino]-5-(trifluoromethyl)pyrimidin-4-yl]-1~{H}-indole-6-carbonitrile"
    common_name = "SY-5609"
    target_protein = "Cyclin-dependent kinase 7 (CDK7)"
    
    # The IC50 value from literature (PMID: 34559569) is used as the binding affinity.
    ic50_value_nM = 0.45
    
    print(f"Compound: {compound_name} (also known as {common_name})")
    print(f"Target: {target_protein}")
    print(f"Reported IC50 value: {ic50_value_nM} nM\n")
    
    # Step 2: Define the answer choices and their ranges in nM.
    # Note: 1 uM = 1,000 nM; 1 mM = 1,000,000 nM
    choices = {
        "A": {"range_str": "< 0.1 nM", "min_nM": -float('inf'), "max_nM": 0.1},
        "B": {"range_str": "0.1 - 100 nM", "min_nM": 0.1, "max_nM": 100},
        "C": {"range_str": "0.1 - 100 uM", "min_nM": 100, "max_nM": 100_000},
        "D": {"range_str": "0.1 - 100 mM", "min_nM": 100_000, "max_nM": 100_000_000},
        "E": {"range_str": "> 100 mM", "min_nM": 100_000_000, "max_nM": float('inf')}
    }
    
    # Step 3: Iterate through choices to find the correct range.
    print("Checking which range the value falls into...")
    correct_choice = None
    for choice, details in choices.items():
        if details["min_nM"] <= ic50_value_nM < details["max_nM"]:
            correct_choice = choice
            print(f"The value {ic50_value_nM} nM is within the range for Choice {choice}: {details['range_str']}.")
            
            # Step 4: Display the final equation as requested.
            print("\nThe final equation showing this relationship is:")
            # We output each number in the final equation
            print(f"{details['min_nM']} nM <= {ic50_value_nM} nM < {details['max_nM']} nM")
            break
            
    if not correct_choice:
        print("Could not find a matching range for the given value.")

# Execute the function
find_binding_affinity_range()