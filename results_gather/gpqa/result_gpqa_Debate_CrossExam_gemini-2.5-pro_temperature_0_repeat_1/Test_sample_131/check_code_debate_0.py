import math
from collections import Counter

def check_nmr_answer():
    """
    Checks the correctness of the proposed answer by simulating the NMR spectrum
    of the mixed compounds and comparing it to the experimental data.
    """
    # Predicted 1H NMR data for each compound based on chemical principles.
    # A key detail is the multiplicity of 1,2,3,4-tetramethylbenzene's aromatic protons.
    compounds_data = {
        "1,2,4,5-tetramethylbenzene": {
            "aromatic": {"integrations": [2], "multiplicities": ["singlet"]},
            "alkyl": {"integrations": [12], "multiplicities": ["singlet"]}
        },
        "1,2,3,5-tetramethylbenzene": {
            "aromatic": {"integrations": [1, 1], "multiplicities": ["singlet", "singlet"]},
            "alkyl": {"integrations": [6, 3, 3], "multiplicities": ["singlet", "singlet", "singlet"]}
        },
        "1,2,3,4-tetramethylbenzene": {
            "aromatic": {"integrations": [2], "multiplicities": ["multiplet"]},
            "alkyl": {"integrations": [6, 6], "multiplicities": ["singlet", "singlet"]}
        },
        "1,4-diethylbenzene": {
            "aromatic": {"integrations": [4], "multiplicities": ["singlet"]},
            "alkyl": {"integrations": [4, 6], "multiplicities": ["quartet", "triplet"]}
        }
    }

    # Options provided in the question
    options = {
        "A": ["1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "B": ["1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "C": ["1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"],
        "D": ["1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"]
    }

    # Experimental data from the question
    experimental_data = {
        "aromatic_signals": 2,
        "aromatic_types": ["singlet", "singlet"],
        "aromatic_ratio": [1, 1],
        "alkyl_signals": 3,
        "alkyl_types": ["singlet", "singlet", "singlet"],
        "alkyl_ratio": [2, 1, 1]
    }

    # The answer to check
    answer_to_check = "D"
    
    # --- Helper function to simplify integration ratios ---
    def simplify_ratio(integrations):
        if not integrations:
            return []
        common_divisor = integrations[0]
        for i in range(1, len(integrations)):
            common_divisor = math.gcd(common_divisor, integrations[i])
        return [i // common_divisor for i in integrations]

    # --- Main checking logic ---
    compound1_name, compound2_name = options[answer_to_check]
    compound1 = compounds_data[compound1_name]
    compound2 = compounds_data[compound2_name]

    # Combine predicted signals for a 1:1 mixture
    combined_alkyl_integrations = compound1["alkyl"]["integrations"] + compound2["alkyl"]["integrations"]
    combined_alkyl_multiplicities = compound1["alkyl"]["multiplicities"] + compound2["alkyl"]["multiplicities"]
    
    combined_aromatic_integrations = compound1["aromatic"]["integrations"] + compound2["aromatic"]["integrations"]
    combined_aromatic_multiplicities = compound1["aromatic"]["multiplicities"] + compound2["aromatic"]["multiplicities"]

    # Check Alkyl Region
    if len(combined_alkyl_integrations) != experimental_data["alkyl_signals"]:
        return f"Incorrect. The mixture produces {len(combined_alkyl_integrations)} alkyl signals, but the spectrum shows {experimental_data['alkyl_signals']}."
    if Counter(combined_alkyl_multiplicities) != Counter(experimental_data["alkyl_types"]):
        return f"Incorrect. The mixture produces non-singlet signals in the alkyl region ({combined_alkyl_multiplicities}), but the spectrum shows only singlets."
    if sorted(simplify_ratio(combined_alkyl_integrations)) != sorted(experimental_data["alkyl_ratio"]):
        return f"Incorrect. The mixture produces an alkyl integration ratio of {simplify_ratio(combined_alkyl_integrations)}, not {experimental_data['alkyl_ratio']}."

    # Check Aromatic Region
    if len(combined_aromatic_integrations) != experimental_data["aromatic_signals"]:
        return f"Incorrect. The mixture produces {len(combined_aromatic_integrations)} aromatic signals, but the spectrum shows {experimental_data['aromatic_signals']}."
    if sorted(simplify_ratio(combined_aromatic_integrations)) != sorted(experimental_data["aromatic_ratio"]):
        return f"Incorrect. The mixture produces an aromatic integration ratio of {simplify_ratio(combined_aromatic_integrations)}, not {experimental_data['aromatic_ratio']}."
    if Counter(combined_aromatic_multiplicities) != Counter(experimental_data["aromatic_types"]):
        return (f"Incorrect. The answer does not satisfy the constraint of having two singlets in the aromatic region. "
                f"The predicted mixture of {compound1_name} and {compound2_name} would show signals with multiplicities of {combined_aromatic_multiplicities}. "
                f"Specifically, the two aromatic protons of 1,2,3,4-tetramethylbenzene are adjacent and would couple to form a multiplet, not a singlet.")

    return "Correct"

# Execute the check and print the result
result = check_nmr_answer()
print(result)