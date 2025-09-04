import re

def check_answer():
    """
    Checks the correctness of the final answer for the given organic chemistry problem.
    """
    # The final answer provided by the LLM to be checked.
    final_answer = "<<<B>>>"
    
    # --- Problem Constraints ---
    # 1. Carbon Count Analysis
    # Starting material: 3,4-dimethylhexanedial (6 from hexane + 2 from dimethyl = 8 carbons)
    # Reagent 2: CH3CH2MgBr adds an ethyl group (+2 carbons)
    # Expected final carbon count = 8 + 2 = 10
    expected_carbon_count = 10

    # 2. Functional Group Analysis
    # Reagent 4: O3, H2O is an oxidative ozonolysis.
    # A C=C-H bond is cleaved to a carboxylic acid (-COOH).
    # The highest priority functional group should be a carboxylic acid.
    expected_primary_group = "Carboxylic Acid"

    # --- Options from the question ---
    options = {
        "A": "4,5-dimethylnonane-2,6,7-trione",
        "B": "3,4-dimethyl-5,6-dioxooctanoic acid",
        "C": "4,5-dimethylnonane-2,6,7-trione", # Duplicate of A in the prompt
        "D": "3,4-dimethyl-5,6-dioxooctanal"
    }

    # --- Parsing the LLM's Answer ---
    match = re.search(r'<<<([A-D])>>>', final_answer)
    if not match:
        return "Incorrect format: The final answer is not in the format <<<A>>>, <<<B>>>, etc."
    
    chosen_option_letter = match.group(1)
    chosen_name = options.get(chosen_option_letter)

    if not chosen_name:
        return f"Invalid option: The chosen option '{chosen_option_letter}' is not one of the valid choices A, B, C, or D."

    # --- Helper functions to parse IUPAC names ---
    def get_carbon_count(name):
        count = 0
        # Parent chain
        if "octan" in name: count += 8
        elif "nonan" in name: count += 9
        # Substituents
        if "dimethyl" in name: count += 2
        elif "methyl" in name: count += 1
        return count

    def get_primary_functional_group(name):
        if name.endswith("oic acid"):
            return "Carboxylic Acid"
        if name.endswith("al"):
            return "Aldehyde"
        if name.endswith("one"):
            return "Ketone"
        return "Unknown"

    # --- Checking the chosen answer against constraints ---
    actual_carbon_count = get_carbon_count(chosen_name)
    actual_primary_group = get_primary_functional_group(chosen_name)

    # Check Constraint 1: Carbon Count
    if actual_carbon_count != expected_carbon_count:
        return (f"Incorrect. The carbon count is wrong. "
                f"The starting material has 8 carbons and the Grignard reaction adds 2, "
                f"so the final product must have {expected_carbon_count} carbons. "
                f"The chosen answer '{chosen_name}' has {actual_carbon_count} carbons.")

    # Check Constraint 2: Functional Group
    if actual_primary_group != expected_primary_group:
        return (f"Incorrect. The primary functional group is wrong. "
                f"The final step is an oxidative ozonolysis (O3, H2O), which must produce a {expected_primary_group}. "
                f"The chosen answer '{chosen_name}' is an {actual_primary_group}.")

    # If all constraints are satisfied
    return "Correct"

# Run the check
result = check_answer()
print(result)