import math
from functools import reduce

def simplify_ratio(numbers):
    """Simplifies a list of numbers to their simplest integer ratio."""
    if not numbers or len(numbers) == 0:
        return []
    # Find the greatest common divisor of all numbers in the list
    try:
        divisor = reduce(math.gcd, numbers)
    except TypeError:
        # Handle case where numbers might not be integers
        return None
    # Divide each number by the GCD
    return [n // divisor for n in numbers]

def check_correctness():
    """
    Checks the correctness of the LLM's answer by simulating the 1H NMR spectrum
    of the proposed mixture and comparing it to the experimental data.
    """
    # Define the NMR properties of each individual compound based on chemical principles
    compounds_db = {
        "1,2,4,5-tetramethylbenzene": {
            "aromatic": {"signals": 1, "types": ["singlet"], "integrations": [2]},
            "alkyl": {"signals": 1, "types": ["singlet"], "integrations": [12]}
        },
        "1,2,3,5-tetramethylbenzene": {
            "aromatic": {"signals": 2, "types": ["singlet", "singlet"], "integrations": [1, 1]},
            "alkyl": {"signals": 3, "types": ["singlet", "singlet", "singlet"], "integrations": [6, 3, 3]}
        },
        "1,2,3,4-tetramethylbenzene": {
            "aromatic": {"signals": 1, "types": ["singlet"], "integrations": [2]},
            "alkyl": {"signals": 2, "types": ["singlet", "singlet"], "integrations": [6, 6]}
        },
        "1,4-diethylbenzene": {
            "aromatic": {"signals": 1, "types": ["singlet"], "integrations": [4]},
            "alkyl": {"signals": 2, "types": ["quartet", "triplet"], "integrations": [4, 6]}
        }
    }

    # Define the options from the question
    options = {
        "A": ["1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"],
        "B": ["1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"],
        "C": ["1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "D": ["1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"]
    }

    # The final answer provided by the LLM to be checked
    llm_answer = "A"

    # The target spectrum from the question
    target_spectrum = {
        "aromatic": {"signals": 2, "ratio": [1, 1]},
        "alkyl": {"signals": 3, "ratio": [2, 1, 1]}
    }

    # --- Verification Logic ---
    try:
        compound_names = options[llm_answer]
        c1_name, c2_name = compound_names[0], compound_names[1]
        c1_data = compounds_db[c1_name]
        c2_data = compounds_db[c2_name]
    except KeyError:
        return f"Invalid answer option '{llm_answer}'. Valid options are A, B, C, D."

    # 1. Check if all alkyl signals are singlets, a key constraint
    combined_alkyl_types = c1_data["alkyl"]["types"] + c2_data["alkyl"]["types"]
    if not all(t == "singlet" for t in combined_alkyl_types):
        non_singlet_compound = c1_name if "singlet" not in c1_data["alkyl"]["types"] else c2_name
        return f"Incorrect. The mixture contains {non_singlet_compound}, which does not produce only singlets in the alkyl region. The problem states all alkyl signals are singlets."

    # 2. Check Aromatic Region
    combined_aromatic_signals = c1_data["aromatic"]["signals"] + c2_data["aromatic"]["signals"]
    if combined_aromatic_signals != target_spectrum["aromatic"]["signals"]:
        return f"Incorrect. The mixture produces {combined_aromatic_signals} aromatic signals, but the spectrum shows {target_spectrum['aromatic']['signals']}."

    combined_aromatic_integrations = c1_data["aromatic"]["integrations"] + c2_data["aromatic"]["integrations"]
    calculated_aromatic_ratio = simplify_ratio(combined_aromatic_integrations)
    if sorted(calculated_aromatic_ratio) != sorted(target_spectrum["aromatic"]["ratio"]):
        return f"Incorrect. The aromatic signal ratio is {calculated_aromatic_ratio} (simplified), which does not match the target ratio of {target_spectrum['aromatic']['ratio']}."

    # 3. Check Alkyl Region
    combined_alkyl_signals = c1_data["alkyl"]["signals"] + c2_data["alkyl"]["signals"]
    if combined_alkyl_signals != target_spectrum["alkyl"]["signals"]:
        return f"Incorrect. The mixture produces {combined_alkyl_signals} alkyl signals, but the spectrum shows {target_spectrum['alkyl']['signals']}."

    combined_alkyl_integrations = c1_data["alkyl"]["integrations"] + c2_data["alkyl"]["integrations"]
    calculated_alkyl_ratio = simplify_ratio(combined_alkyl_integrations)
    if sorted(calculated_alkyl_ratio) != sorted(target_spectrum["alkyl"]["ratio"]):
        return f"Incorrect. The alkyl signal ratio is {calculated_alkyl_ratio} (simplified), which does not match the target ratio of {target_spectrum['alkyl']['ratio']}."

    # If all checks pass, the answer is correct
    return "Correct"

# The code block above can be executed to check the answer.
# print(check_correctness())