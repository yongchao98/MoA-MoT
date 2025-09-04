import math
from collections import Counter

def get_simplified_ratio(integrations):
    """
    Calculates the simplest integer ratio for a list of integrations.
    For example, [12, 6, 6] becomes [2, 1, 1].
    """
    if not integrations:
        return []
    
    # Find the greatest common divisor (GCD) of all integrations
    common_divisor = integrations[0]
    for i in range(1, len(integrations)):
        common_divisor = math.gcd(common_divisor, integrations[i])
    
    # Divide each integration by the GCD to get the simplest ratio
    ratio = [i // common_divisor for i in integrations]
    
    # Return the ratio sorted in descending order for consistent comparison
    return sorted(ratio, reverse=True)

def check_correctness():
    """
    Checks the correctness of the LLM's answer by modeling the NMR spectra.
    """
    # Step 1: Define the predicted spectra for each individual compound
    compounds_db = {
        "1,2,4,5-tetramethylbenzene": {
            "aromatic": {"singlet": [2]},
            "alkyl": {"singlet": [12]}
        },
        "1,2,3,5-tetramethylbenzene": {
            "aromatic": {"singlet": [1, 1]},
            "alkyl": {"singlet": [6, 3, 3]}
        },
        "1,2,3,4-tetramethylbenzene": {
            "aromatic": {"singlet": [2]},
            "alkyl": {"singlet": [6, 6]}
        },
        "1,4-diethylbenzene": {
            "aromatic": {"singlet": [4]},
            "alkyl": {"quartet": [4], "triplet": [6]}
        }
    }

    # Step 2: Define the options from the question
    options = {
        "A": ["1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"],
        "B": ["1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "C": ["1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"],
        "D": ["1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"]
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer = "A"

    # Step 3: Define the target spectrum from the question's description
    target_spectrum = {
        "aromatic": {
            "signals": 2,
            "multiplicity": "singlet",
            "ratio": [1, 1]
        },
        "alkyl": {
            "signals": 3,
            "multiplicity": "singlet",
            "ratio": [2, 1, 1]
        }
    }

    # Step 4: Iterate through each option and check if it matches the target spectrum
    matching_option = None
    for option_letter, compound_names in options.items():
        c1_name, c2_name = compound_names
        c1_data = compounds_db[c1_name]
        c2_data = compounds_db[c2_name]

        # --- Check Alkyl Region Constraints ---
        # Constraint: All signals must be singlets
        if target_spectrum['alkyl']['multiplicity'] not in c1_data['alkyl'] or \
           target_spectrum['alkyl']['multiplicity'] not in c2_data['alkyl']:
            continue

        # Combine alkyl signals from both compounds
        alkyl_integrations = c1_data['alkyl']['singlet'] + c2_data['alkyl']['singlet']
        
        # Constraint: Check number of signals
        if len(alkyl_integrations) != target_spectrum['alkyl']['signals']:
            continue

        # Constraint: Check integration ratio
        alkyl_ratio = get_simplified_ratio(alkyl_integrations)
        if alkyl_ratio != sorted(target_spectrum['alkyl']['ratio'], reverse=True):
            continue

        # --- Check Aromatic Region Constraints ---
        # Combine aromatic signals
        aromatic_integrations = c1_data['aromatic']['singlet'] + c2_data['aromatic']['singlet']

        # Constraint: Check number of signals
        if len(aromatic_integrations) != target_spectrum['aromatic']['signals']:
            continue

        # Constraint: Check integration ratio
        aromatic_ratio = get_simplified_ratio(aromatic_integrations)
        if aromatic_ratio != sorted(target_spectrum['aromatic']['ratio'], reverse=True):
            continue

        # If all checks pass, this is the correct option
        matching_option = option_letter
        break

    # Step 5: Compare the identified correct option with the LLM's answer
    if matching_option is None:
        return "Analysis failed: No option perfectly matches the spectral data based on a straightforward interpretation. The problem may require assuming accidental signal overlaps, which the provided answer did not state."

    if llm_answer == matching_option:
        return "Correct"
    else:
        return f"Incorrect. The provided answer is {llm_answer}, but the correct option is {matching_option}. The mixture of {options[matching_option][0]} and {options[matching_option][1]} is the only one that produces {target_spectrum['aromatic']['signals']} aromatic singlets with a 1:1 ratio and {target_spectrum['alkyl']['signals']} alkyl singlets with a 2:1:1 ratio without requiring assumptions of accidental signal overlap."

# Execute the check
result = check_correctness()
print(result)