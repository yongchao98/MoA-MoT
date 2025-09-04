import math
from functools import reduce
from collections import Counter

def check_answer():
    """
    Checks the correctness of the answer to the NMR spectroscopy question.
    """

    # --- Step 1: Define the spectral properties of each candidate compound ---
    # Each compound is represented as a dictionary containing its formula and
    # predicted 1H NMR signals. Signals are tuples of (integration, multiplicity).
    # s=singlet, q=quartet, t=triplet
    compounds = {
        "1,2,4,5-tetramethylbenzene": {
            "formula": "C10H14",
            # Aromatic: 2 equivalent protons (H3, H6) -> 1 signal
            # Aliphatic: 4 equivalent methyl groups -> 1 signal
            "aromatic_signals": [(2, 's')],
            "aliphatic_signals": [(12, 's')]
        },
        "1,2,3,5-tetramethylbenzene": {
            "formula": "C10H14",
            # Aromatic: 2 equivalent protons (H4, H6) -> 1 signal
            # Aliphatic: 3 sets of methyl groups: (C1,C3), C2, C5 -> 3 signals
            "aromatic_signals": [(2, 's')],
            "aliphatic_signals": [(6, 's'), (3, 's'), (3, 's')]
        },
        "1,2,3,4-tetramethylbenzene": {
            "formula": "C10H14",
            # Aromatic: 2 equivalent protons (H5, H6) -> 1 signal
            # Aliphatic: 2 sets of methyl groups: (C1,C4), (C2,C3) -> 2 signals
            "aromatic_signals": [(2, 's')],
            "aliphatic_signals": [(6, 's'), (6, 's')]
        },
        "1,4-diethylbenzene": {
            "formula": "C10H14",
            # Aromatic: 4 equivalent protons -> 1 signal
            # Aliphatic: CH2 quartet, CH3 triplet -> 2 signals (not singlets)
            "aromatic_signals": [(4, 's')],
            "aliphatic_signals": [(4, 'q'), (6, 't')]
        }
    }

    # --- Step 2: Define the options and the provided answer ---
    options = {
        "A": ["1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"],
        "B": ["1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"],
        "C": ["1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "D": ["1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"]
    }
    
    llm_answer_key = "B"
    
    # --- Step 3: Define a function to check a given mixture ---
    def check_mixture(compound_names):
        c1_name, c2_name = compound_names
        c1 = compounds[c1_name]
        c2 = compounds[c2_name]

        # Constraint: Molecular Formula
        if not (c1["formula"] == "C10H14" and c2["formula"] == "C10H14"):
            return f"Incorrect formula for {c1_name} or {c2_name}."

        # Combine signals for a 1:1 mixture
        combined_aromatic = c1["aromatic_signals"] + c2["aromatic_signals"]
        combined_aliphatic = c1["aliphatic_signals"] + c2["aliphatic_signals"]

        # Constraint: Aromatic signals (2 singlets, 1:1 ratio)
        if len(combined_aromatic) != 2:
            return f"Aromatic signals check failed: Expected 2 signals, but the mixture of {c1_name} and {c2_name} gives {len(combined_aromatic)}."
        if not all(sig[1] == 's' for sig in combined_aromatic):
            return f"Aromatic signals check failed: Not all signals are singlets for the mixture of {c1_name} and {c2_name}."
        # For a 1:1 mixture, a 1:1 ratio means the integrations must be equal
        if combined_aromatic[0][0] != combined_aromatic[1][0]:
            return f"Aromatic signals check failed: Ratio is not 1:1. Integrations are {combined_aromatic[0][0]}H and {combined_aromatic[1][0]}H."

        # Constraint: Aliphatic signals (3 singlets, 2:1:1 ratio)
        if not all(sig[1] == 's' for sig in combined_aliphatic):
            return f"Aliphatic signals check failed: Not all signals are singlets. The mixture of {c1_name} and {c2_name} contains non-singlet signals."
        if len(combined_aliphatic) != 3:
            return f"Aliphatic signals check failed: Expected 3 signals, but the mixture of {c1_name} and {c2_name} gives {len(combined_aliphatic)}."
        
        # Check ratio
        integrations = sorted([sig[0] for sig in combined_aliphatic], reverse=True)
        if len(integrations) == 0: # Should not happen if count is 3
             return "Aliphatic signals check failed: No integrations found."
        common_divisor = reduce(math.gcd, integrations)
        ratio = [i // common_divisor for i in integrations]
        
        # The target ratio is 2:1:1. We check against a sorted version.
        target_ratio = sorted([2, 1, 1], reverse=True)
        if ratio != target_ratio:
            return f"Aliphatic signals check failed: Ratio is {':'.join(map(str, ratio))}, not 2:1:1."

        return "Correct"

    # --- Step 4: Evaluate the LLM's answer ---
    compounds_to_check = options.get(llm_answer_key)
    if not compounds_to_check:
        return f"Invalid answer key '{llm_answer_key}' provided."
        
    result = check_mixture(compounds_to_check)
    
    return result

# Run the check
correctness_result = check_answer()
print(correctness_result)