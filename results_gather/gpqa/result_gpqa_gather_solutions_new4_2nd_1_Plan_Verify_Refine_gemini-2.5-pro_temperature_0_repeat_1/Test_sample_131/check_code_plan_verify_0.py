import math
from collections import Counter

def check_correctness():
    """
    Checks the correctness of the answer for the given NMR spectroscopy problem.
    """
    # Predicted 1H NMR spectral data for each compound based on symmetry
    spectra = {
        "1,2,4,5-tetramethylbenzene": {
            "aromatic": {"multiplicity": ["singlet"], "integrations": [2]},
            "aliphatic": {"multiplicity": ["singlet"], "integrations": [12]}
        },
        "1,2,3,4-tetramethylbenzene": {
            "aromatic": {"multiplicity": ["singlet"], "integrations": [2]},
            "aliphatic": {"multiplicity": ["singlet", "singlet"], "integrations": [6, 6]}
        },
        "1,2,3,5-tetramethylbenzene": {
            "aromatic": {"multiplicity": ["singlet"], "integrations": [2]},
            "aliphatic": {"multiplicity": ["singlet", "singlet", "singlet"], "integrations": [6, 3, 3]}
        },
        "1,4-diethylbenzene": {
            "aromatic": {"multiplicity": ["singlet"], "integrations": [4]},
            "aliphatic": {"multiplicity": ["quartet", "triplet"], "integrations": [4, 6]}
        }
    }

    # Options provided in the question
    options = {
        "A": ["1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"],
        "B": ["1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "C": ["1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "D": ["1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"]
    }

    # The final answer provided by the LLM to be checked
    llm_answer_key = "A"

    # Target spectrum from the problem description
    target_aromatic_signals = 2
    target_aromatic_ratio = sorted([1, 1])
    target_aliphatic_signals = 3
    target_aliphatic_ratio = sorted([2, 1, 1])

    # Get the compounds from the selected answer
    compounds_to_check = options[llm_answer_key]
    compound1_name, compound2_name = compounds_to_check[0], compounds_to_check[1]
    spec1, spec2 = spectra[compound1_name], spectra[compound2_name]

    # --- Check Aliphatic Region ---
    aliphatic_multiplicity = spec1["aliphatic"]["multiplicity"] + spec2["aliphatic"]["multiplicity"]
    if not all(m == "singlet" for m in aliphatic_multiplicity):
        non_singlet_compound = compound1_name if "singlet" not in spec1["aliphatic"]["multiplicity"] else compound2_name
        return f"Incorrect. The aliphatic signals must all be singlets, but compound '{non_singlet_compound}' produces non-singlet signals (quartets/triplets)."

    aliphatic_integrations = spec1["aliphatic"]["integrations"] + spec2["aliphatic"]["integrations"]
    if len(aliphatic_integrations) != target_aliphatic_signals:
        return f"Incorrect. The mixture produces {len(aliphatic_integrations)} aliphatic signals, but the question specifies {target_aliphatic_signals} signals."

    # Check ratio by normalizing with the smallest integration unit
    integrations_copy = sorted(aliphatic_integrations)
    unit = integrations_copy[0] # Smallest integration value
    observed_ratio = sorted([round(i / unit) for i in integrations_copy])
    if observed_ratio != target_aliphatic_ratio:
        return f"Incorrect. The aliphatic signal integration ratio is {':'.join(map(str, sorted(aliphatic_integrations, reverse=True)))}, which simplifies to {':'.join(map(str, observed_ratio))}, not the required 2:1:1 ratio."

    # --- Check Aromatic Region ---
    aromatic_multiplicity = spec1["aromatic"]["multiplicity"] + spec2["aromatic"]["multiplicity"]
    if not all(m == "singlet" for m in aromatic_multiplicity):
        return "Incorrect. The aromatic signals must all be singlets."

    aromatic_integrations = spec1["aromatic"]["integrations"] + spec2["aromatic"]["integrations"]
    if len(aromatic_integrations) != target_aromatic_signals:
        return f"Incorrect. The mixture produces {len(aromatic_integrations)} aromatic signals, but the question specifies {target_aromatic_signals} signals."

    # Check ratio
    integrations_copy = sorted(aromatic_integrations)
    unit = integrations_copy[0]
    observed_ratio = sorted([round(i / unit) for i in integrations_copy])
    if observed_ratio != target_aromatic_ratio:
        return f"Incorrect. The aromatic signal integration ratio is not the required 1:1 ratio."

    return "Correct"

# Run the check
result = check_correctness()
print(result)