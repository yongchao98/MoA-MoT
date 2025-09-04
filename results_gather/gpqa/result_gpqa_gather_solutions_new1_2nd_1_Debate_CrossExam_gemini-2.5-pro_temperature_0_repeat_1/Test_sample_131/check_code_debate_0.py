import math
from collections import Counter

def check_answer():
    """
    Checks the correctness of the LLM's answer by simulating the NMR spectrum of the proposed mixture.
    """

    # Step 1: Define the predicted 1H NMR data for each individual compound.
    # The data is stored as a dictionary where each signal is a tuple of (integration, multiplicity).
    # Integrations are based on the number of protons (H).
    nmr_data = {
        "1,2,4,5-tetramethylbenzene": {
            "aromatic": [(2, "singlet")],  # 2 equivalent aromatic H
            "alkyl": [(12, "singlet")]      # 4 equivalent methyl groups (4*3=12H)
        },
        "1,2,3,5-tetramethylbenzene": {
            "aromatic": [(1, "singlet"), (1, "singlet")], # 2 non-equivalent aromatic H
            "alkyl": [(6, "singlet"), (3, "singlet"), (3, "singlet")] # 3 sets of methyl groups (2, 1, 1)
        },
        "1,2,3,4-tetramethylbenzene": {
            "aromatic": [(2, "singlet")],  # 2 equivalent aromatic H
            "alkyl": [(6, "singlet"), (6, "singlet")] # 2 sets of methyl groups (2 each)
        },
        "1,4-diethylbenzene": {
            "aromatic": [(4, "singlet")],  # 4 equivalent aromatic H
            "alkyl": [(4, "quartet"), (6, "triplet")] # 2xCH2 (4H), 2xCH3 (6H)
        }
    }

    # Step 2: Define the options from the question.
    options = {
        "A": ["1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"],
        "B": ["1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"],
        "C": ["1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "D": ["1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"]
    }

    # The final answer provided by the LLM.
    llm_answer = "A"

    # Helper function to simplify an integration ratio (e.g., [12, 6, 6] -> [2, 1, 1])
    def simplify_ratio(numbers):
        if not numbers:
            return []
        common_divisor = numbers[0]
        for num in numbers[1:]:
            common_divisor = math.gcd(common_divisor, num)
        return [num // common_divisor for num in numbers]

    # Step 3: Get the compounds for the proposed answer and check them against the constraints.
    compounds_to_check = options.get(llm_answer)

    if compounds_to_check is None:
        return f"Invalid answer option '{llm_answer}' provided."

    # Combine the signals from the two compounds in the proposed mixture.
    combined_aromatic_signals = []
    combined_alkyl_signals = []
    for c_name in compounds_to_check:
        data = nmr_data[c_name]
        combined_aromatic_signals.extend(data["aromatic"])
        combined_alkyl_signals.extend(data["alkyl"])

    # --- Check Aromatic Region Constraints ---
    # Constraint: Two signals
    if len(combined_aromatic_signals) != 2:
        return f"Incorrect. The proposed mixture (Option {llm_answer}) produces {len(combined_aromatic_signals)} aromatic signals, but the spectrum requires 2."
    
    # Constraint: All singlets
    for _, multiplicity in combined_aromatic_signals:
        if multiplicity != "singlet":
            return f"Incorrect. The proposed mixture (Option {llm_answer}) produces a {multiplicity} in the aromatic region, but all signals must be singlets."

    # Constraint: 1:1 ratio
    aromatic_integrations = sorted([s[0] for s in combined_aromatic_signals], reverse=True)
    simplified_aromatic_ratio = simplify_ratio(aromatic_integrations)
    if simplified_aromatic_ratio != [1, 1]:
        return f"Incorrect. The aromatic signal integration ratio for Option {llm_answer} is {':'.join(map(str, simplified_aromatic_ratio))}, not 1:1."

    # --- Check Alkyl Region Constraints ---
    # Constraint: All singlets
    for _, multiplicity in combined_alkyl_signals:
        if multiplicity != "singlet":
            return f"Incorrect. The proposed mixture (Option {llm_answer}) contains a compound that produces a {multiplicity} in the alkyl region, but all signals must be singlets."

    # Constraint: Three signals
    if len(combined_alkyl_signals) != 3:
        return f"Incorrect. The proposed mixture (Option {llm_answer}) produces {len(combined_alkyl_signals)} alkyl signals, but the spectrum requires 3."

    # Constraint: 2:1:1 ratio
    alkyl_integrations = sorted([s[0] for s in combined_alkyl_signals], reverse=True)
    simplified_alkyl_ratio = simplify_ratio(alkyl_integrations)
    if simplified_alkyl_ratio != [2, 1, 1]:
        return f"Incorrect. The alkyl signal integration ratio for Option {llm_answer} is {':'.join(map(str, simplified_alkyl_ratio))}, not 2:1:1."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
print(check_answer())