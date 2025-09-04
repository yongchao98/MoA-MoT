import collections

def check_answer():
    """
    This function checks the correctness of the provided answer to the NMR spectroscopy question.
    It simulates the 1H NMR spectrum for each possible mixture and compares it to the experimental data.
    """

    # 1. Define the predicted 1H NMR properties for each individual compound.
    # This is based on chemical structure and symmetry.
    # 'integration' refers to the number of protons (H).
    compounds_db = {
        "1,2,4,5-tetramethylbenzene": {
            "aromatic": [{"multiplicity": "singlet", "integration": 2}],
            "aliphatic": [{"multiplicity": "singlet", "integration": 12}]
        },
        "1,2,3,5-tetramethylbenzene": {
            # Aromatic protons at C4 and C6 are not equivalent.
            "aromatic": [{"multiplicity": "singlet", "integration": 1}, {"multiplicity": "singlet", "integration": 1}],
            # Three sets of methyl groups: (C1,C3), C2, C5.
            "aliphatic": [{"multiplicity": "singlet", "integration": 6}, {"multiplicity": "singlet", "integration": 3}, {"multiplicity": "singlet", "integration": 3}]
        },
        "1,2,3,4-tetramethylbenzene": {
            # Aromatic protons at C5 and C6 are equivalent (AA' system), treated as a singlet as per the problem's context.
            "aromatic": [{"multiplicity": "singlet", "integration": 2}],
            # Two sets of methyl groups: (C1,C4) and (C2,C3).
            "aliphatic": [{"multiplicity": "singlet", "integration": 6}, {"multiplicity": "singlet", "integration": 6}]
        },
        "1,4-diethylbenzene": {
            "aromatic": [{"multiplicity": "singlet", "integration": 4}],
            # Ethyl groups show coupling.
            "aliphatic": [{"multiplicity": "quartet", "integration": 4}, {"multiplicity": "triplet", "integration": 6}]
        }
    }

    # 2. Define the experimental data from the question.
    expected_spectrum = {
        "aromatic": {"signals": 2, "multiplicity": "singlet", "ratio": [1, 1]},
        "aliphatic": {"signals": 3, "multiplicity": "singlet", "ratio": [2, 1, 1]}
    }

    # 3. Define the options given in the question.
    options = {
        "A": ["1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"],
        "B": ["1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "C": ["1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"],
        "D": ["1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"]
    }
    
    # The final answer provided by the LLM to be checked.
    llm_answer = "A"

    # Helper function to check if integration values match a given ratio.
    def check_ratio(integrations, ratio):
        if len(integrations) != len(ratio):
            return False
        
        sorted_integrations = sorted(integrations, reverse=True)
        sorted_ratio = sorted(ratio, reverse=True)
        
        # Normalize both lists by their smallest value and compare.
        base_integration = sorted_integrations[-1]
        base_ratio = sorted_ratio[-1]
        
        if base_integration == 0 or base_ratio == 0:
            return False # Avoid division by zero

        normalized_integrations = [i / base_integration for i in sorted_integrations]
        normalized_ratio = [r / base_ratio for r in sorted_ratio]
        
        # Compare the normalized lists with a small tolerance for floating point errors.
        return all(abs(a - b) < 1e-9 for a, b in zip(normalized_integrations, normalized_ratio))

    # 4. Iterate through each option and check if it matches the experimental data.
    correct_option = None
    for option_key, compound_names in options.items():
        comp1_name, comp2_name = compound_names
        comp1 = compounds_db[comp1_name]
        comp2 = compounds_db[comp2_name]

        # --- Check Aliphatic Region ---
        aliphatic_signals = comp1["aliphatic"] + comp2["aliphatic"]
        
        # Constraint: All signals must be singlets.
        if not all(s["multiplicity"] == "singlet" for s in aliphatic_signals):
            continue # This option is incorrect, move to the next one.
            
        # Constraint: Must have the correct number of signals.
        if len(aliphatic_signals) != expected_spectrum["aliphatic"]["signals"]:
            continue
            
        # Constraint: Must have the correct integration ratio.
        aliphatic_integrations = [s["integration"] for s in aliphatic_signals]
        if not check_ratio(aliphatic_integrations, expected_spectrum["aliphatic"]["ratio"]):
            continue

        # --- Check Aromatic Region ---
        aromatic_signals = comp1["aromatic"] + comp2["aromatic"]
        
        # Constraint: All signals must be singlets.
        if not all(s["multiplicity"] == "singlet" for s in aromatic_signals):
            continue
            
        # Constraint: Must have the correct number of signals.
        if len(aromatic_signals) != expected_spectrum["aromatic"]["signals"]:
            continue
            
        # Constraint: Must have the correct integration ratio.
        aromatic_integrations = [s["integration"] for s in aromatic_signals]
        if not check_ratio(aromatic_integrations, expected_spectrum["aromatic"]["ratio"]):
            continue
            
        # If all checks pass, this is the correct option.
        correct_option = option_key
        break # Found the correct option, no need to check others.

    # 5. Compare the identified correct option with the LLM's answer.
    if correct_option == llm_answer:
        return "Correct"
    elif correct_option is None:
        return f"The provided answer '{llm_answer}' is incorrect. No option perfectly matches the experimental data based on the analysis. The analysis for '{llm_answer}' failed because its predicted spectrum does not match the given constraints."
    else:
        return f"The provided answer '{llm_answer}' is incorrect. The correct answer should be '{correct_option}'. Option '{llm_answer}' fails the check, while option '{correct_option}' perfectly matches all spectral data (aromatic: 2 singlets in 1:1 ratio; aliphatic: 3 singlets in 2:1:1 ratio)."

# Execute the check and print the result.
print(check_answer())