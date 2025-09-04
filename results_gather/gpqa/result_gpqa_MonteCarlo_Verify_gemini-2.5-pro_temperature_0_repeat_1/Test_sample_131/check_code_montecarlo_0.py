import math
from collections import Counter

def get_simplified_ratio(integrals):
    """Calculates the simplest integer ratio for a list of numbers."""
    if not integrals:
        return []
    # Find the greatest common divisor (GCD) of all integrals
    common_divisor = integrals[0]
    for i in range(1, len(integrals)):
        common_divisor = math.gcd(common_divisor, integrals[i])
    
    # Divide each integral by the GCD and sort for consistent comparison
    return sorted([i // common_divisor for i in integrals])

def check_answer():
    """
    Checks if the proposed mixture of compounds matches the NMR data from the question.
    """
    # --- Ground Truth Data ---
    # Define the predicted 1H NMR singlet signals for each compound.
    # Format: { 'aromatic': [integrations], 'aliphatic': [integrations] }
    # Non-singlet signals are ignored as per the question's focus.
    compound_spectra = {
        "1,2,4,5-tetramethylbenzene": {
            "aromatic": [2],
            "aliphatic": [12]
        },
        "1,2,3,5-tetramethylbenzene": {
            "aromatic": [2],
            "aliphatic": [6, 3, 3]
        },
        "1,2,3,4-tetramethylbenzene": {
            "aromatic": [2],
            "aliphatic": [6, 6]
        },
        "1,4-diethylbenzene": {
            "aromatic": [4],
            "aliphatic": []  # Has non-singlet signals (quartet, triplet)
        }
    }

    # --- Problem Constraints ---
    # Define the target spectrum based on the question.
    target_spectrum = {
        "aromatic_singlet_count": 2,
        "aromatic_ratio": [1, 1],
        "aliphatic_singlet_count": 3,
        "aliphatic_ratio": [1, 1, 2]  # Sorted version of 2:1:1
    }

    # --- Proposed Answer ---
    # This corresponds to option B from the multiple-choice question.
    proposed_answer_compounds = ["1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"]
    
    # --- Verification Logic ---
    c1_name, c2_name = proposed_answer_compounds
    c1_data = compound_spectra[c1_name]
    c2_data = compound_spectra[c2_name]

    # 1. Check Aromatic Signals
    combined_aromatic_integrals = c1_data["aromatic"] + c2_data["aromatic"]
    
    if len(combined_aromatic_integrals) != target_spectrum["aromatic_singlet_count"]:
        return (f"Incorrect. The proposed mixture produces {len(combined_aromatic_integrals)} aromatic singlets, "
                f"but the question specifies {target_spectrum['aromatic_singlet_count']}.")

    aromatic_ratio = get_simplified_ratio(combined_aromatic_integrals)
    if aromatic_ratio != target_spectrum["aromatic_ratio"]:
        return (f"Incorrect. The aromatic signal ratio is {aromatic_ratio}, "
                f"but the question specifies a 1:1 ratio ({target_spectrum['aromatic_ratio']}).")

    # 2. Check Aliphatic Signals
    combined_aliphatic_integrals = c1_data["aliphatic"] + c2_data["aliphatic"]

    if len(combined_aliphatic_integrals) != target_spectrum["aliphatic_singlet_count"]:
        return (f"Incorrect. The proposed mixture produces {len(combined_aliphatic_integrals)} aliphatic singlets, "
                f"but the question specifies {target_spectrum['aliphatic_singlet_count']}.")

    aliphatic_ratio = get_simplified_ratio(combined_aliphatic_integrals)
    if aliphatic_ratio != target_spectrum["aliphatic_ratio"]:
        return (f"Incorrect. The aliphatic signal ratio is {aliphatic_ratio} (from integrals {combined_aliphatic_integrals}), "
                f"but the question specifies a 2:1:1 ratio ({target_spectrum['aliphatic_ratio']}).")

    # 3. If all checks pass
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)