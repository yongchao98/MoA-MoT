import math
from collections import Counter

def check_correctness():
    """
    This function checks the correctness of the LLM's answer by simulating the 1H NMR spectrum for each option.
    It verifies the number of signals, their multiplicities, and their integration ratios against the problem's constraints.
    """

    # Step 1: Define the predicted 1H NMR data for each compound based on its structure and symmetry.
    # Data is stored as: { "aromatic": [(multiplicity, integration), ...], "aliphatic": [...] }
    nmr_predictions = {
        "1,2,4,5-tetramethylbenzene": {
            "formula": "C10H14",
            "aromatic": [("singlet", 2)],  # 2 equivalent aromatic protons
            "aliphatic": [("singlet", 12)] # 4 equivalent methyl groups
        },
        "1,2,3,5-tetramethylbenzene": {
            "formula": "C10H14",
            "aromatic": [("singlet", 1), ("singlet", 1)], # 2 non-equivalent aromatic protons
            "aliphatic": [("singlet", 3), ("singlet", 3), ("singlet", 3), ("singlet", 3)] # 4 non-equivalent methyl groups
        },
        "1,2,3,4-tetramethylbenzene": {
            "formula": "C10H14",
            "aromatic": [("singlet", 2)],  # 2 equivalent aromatic protons
            "aliphatic": [("singlet", 6), ("singlet", 6)] # 2 pairs of equivalent methyl groups
        },
        "1,4-diethylbenzene": {
            "formula": "C10H14",
            "aromatic": [("singlet", 4)],  # 4 equivalent aromatic protons
            "aliphatic": [("quartet", 4), ("triplet", 6)] # -CH2- quartet and -CH3 triplet
        }
    }

    # Step 2: Define the experimental data from the question as constraints.
    constraints = {
        "formula": "C10H14",
        "aromatic_signals": 2,
        "aromatic_multiplicity": "singlet",
        "aromatic_ratio": Counter([1, 1]),
        "aliphatic_signals": 3,
        "aliphatic_multiplicity": "singlet",
        "aliphatic_ratio": Counter([2, 1, 1])
    }

    # Step 3: Define the options given in the question.
    options = {
        "A": ("1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"),
        "B": ("1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"),
        "C": ("1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"),
        "D": ("1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene")
    }

    llm_answer_key = "C"
    
    # Helper function to simplify integration ratios
    def get_simplified_ratio(integrations):
        if not integrations:
            return []
        common_divisor = integrations[0]
        for i in integrations[1:]:
            common_divisor = math.gcd(common_divisor, i)
        return Counter([i // common_divisor for i in integrations])

    # Step 4: Analyze the pair of compounds from the LLM's answer.
    compound1_name, compound2_name = options[llm_answer_key]
    
    # Check molecular formula constraint
    if nmr_predictions[compound1_name]["formula"] != constraints["formula"] or \
       nmr_predictions[compound2_name]["formula"] != constraints["formula"]:
        return f"Incorrect. The molecular formula for one or both compounds in option {llm_answer_key} is not {constraints['formula']}."

    # Combine the predicted spectra for the 1:1 mixture
    c1_data = nmr_predictions[compound1_name]
    c2_data = nmr_predictions[compound2_name]
    
    combined_aromatic = c1_data["aromatic"] + c2_data["aromatic"]
    combined_aliphatic = c1_data["aliphatic"] + c2_data["aliphatic"]

    # Step 5: Check the combined spectrum against the constraints.
    # Check Aromatic Signals
    if len(combined_aromatic) != constraints["aromatic_signals"]:
        return f"Incorrect. The answer {llm_answer_key} is wrong because the mixture produces {len(combined_aromatic)} aromatic signals, but the question specifies {constraints['aromatic_signals']}."
    
    for multiplicity, _ in combined_aromatic:
        if multiplicity != constraints["aromatic_multiplicity"]:
            return f"Incorrect. The answer {llm_answer_key} is wrong because it produces an aromatic signal with multiplicity '{multiplicity}', but the question specifies only singlets."

    aromatic_integrations = [integ for _, integ in combined_aromatic]
    if get_simplified_ratio(aromatic_integrations) != constraints["aromatic_ratio"]:
        return f"Incorrect. The answer {llm_answer_key} is wrong because the aromatic signal integration ratio is {':'.join(map(str, sorted(aromatic_integrations)))}, which does not simplify to 1:1."

    # Check Aliphatic Signals
    if len(combined_aliphatic) != constraints["aliphatic_signals"]:
        return f"Incorrect. The answer {llm_answer_key} is wrong because the mixture produces {len(combined_aliphatic)} aliphatic signals, but the question specifies {constraints['aliphatic_signals']}."

    for multiplicity, _ in combined_aliphatic:
        if multiplicity != constraints["aliphatic_multiplicity"]:
            return f"Incorrect. The answer {llm_answer_key} is wrong because it produces an aliphatic signal with multiplicity '{multiplicity}', but the question specifies only singlets."

    aliphatic_integrations = [integ for _, integ in combined_aliphatic]
    if get_simplified_ratio(aliphatic_integrations) != constraints["aliphatic_ratio"]:
        actual_ratio_str = ':'.join(map(str, sorted(get_simplified_ratio(aliphatic_integrations).elements(), reverse=True)))
        expected_ratio_str = ':'.join(map(str, sorted(constraints['aliphatic_ratio'].elements(), reverse=True)))
        return f"Incorrect. The answer {llm_answer_key} is wrong because the aliphatic signal integration ratio is {actual_ratio_str}, but the question specifies a ratio of {expected_ratio_str}."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness()
print(result)