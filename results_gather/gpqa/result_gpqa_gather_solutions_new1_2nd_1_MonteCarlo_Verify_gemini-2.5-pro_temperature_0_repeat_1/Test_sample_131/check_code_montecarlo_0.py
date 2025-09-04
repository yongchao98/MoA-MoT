import math
from functools import reduce

def check_correctness():
    """
    This function checks the correctness of the final answer by simulating the
    NMR spectrum of the proposed mixture and comparing it against the given constraints.
    """

    # --- Step 1: Define the predicted NMR data for each compound ---
    # This data is based on chemical structure and symmetry.
    compounds_data = {
        '1,2,4,5-tetramethylbenzene': {
            'aromatic': [{'protons': 2, 'multiplicity': 'singlet'}],
            'alkyl': [{'protons': 12, 'multiplicity': 'singlet'}]
        },
        '1,2,3,5-tetramethylbenzene': {
            'aromatic': [{'protons': 1, 'multiplicity': 'singlet'}, {'protons': 1, 'multiplicity': 'singlet'}],
            'alkyl': [{'protons': 6, 'multiplicity': 'singlet'}, {'protons': 3, 'multiplicity': 'singlet'}, {'protons': 3, 'multiplicity': 'singlet'}]
        },
        '1,2,3,4-tetramethylbenzene': {
            'aromatic': [{'protons': 2, 'multiplicity': 'singlet'}],
            'alkyl': [{'protons': 6, 'multiplicity': 'singlet'}, {'protons': 6, 'multiplicity': 'singlet'}]
        },
        '1,4-diethylbenzene': {
            'aromatic': [{'protons': 4, 'multiplicity': 'singlet'}],
            'alkyl': [{'protons': 2, 'multiplicity': 'quartet'}, {'protons': 3, 'multiplicity': 'triplet'}]
        }
    }

    # --- Step 2: Define the options and the proposed answer ---
    options = {
        "A": ["1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"],
        "B": ["1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"],
        "C": ["1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "D": ["1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"]
    }
    
    # The final answer from the LLM is 'A'.
    llm_answer = "A"

    # --- Step 3: Combine the spectra for the proposed mixture ---
    compound_names = options[llm_answer]
    c1_data = compounds_data[compound_names[0]]
    c2_data = compounds_data[compound_names[1]]

    combined_aromatic_signals = c1_data['aromatic'] + c2_data['aromatic']
    combined_alkyl_signals = c1_data['alkyl'] + c2_data['alkyl']

    # --- Step 4: Check all constraints from the question ---

    # Aromatic Region Constraints
    if len(combined_aromatic_signals) != 2:
        return f"Incorrect. The mixture for option {llm_answer} has {len(combined_aromatic_signals)} aromatic signals, but the question requires 2."
    if not all(s['multiplicity'] == 'singlet' for s in combined_aromatic_signals):
        return f"Incorrect. The mixture for option {llm_answer} has non-singlet aromatic signals, but the question requires only singlets."
    aromatic_protons = [s['protons'] for s in combined_aromatic_signals]
    if aromatic_protons[0] != aromatic_protons[1]:
        return f"Incorrect. The aromatic signals for option {llm_answer} have a proton ratio of {aromatic_protons[0]}:{aromatic_protons[1]}, but the question requires a 1:1 ratio."

    # Alkyl Region Constraints
    if len(combined_alkyl_signals) != 3:
        return f"Incorrect. The mixture for option {llm_answer} has {len(combined_alkyl_signals)} alkyl signals, but the question requires 3."
    if not all(s['multiplicity'] == 'singlet' for s in combined_alkyl_signals):
        return f"Incorrect. The mixture for option {llm_answer} has non-singlet alkyl signals, but the question requires only singlets."
    
    # Check 2:1:1 ratio
    alkyl_protons = [s['protons'] for s in combined_alkyl_signals]
    def gcd_list(numbers):
        return reduce(math.gcd, numbers)
    common_divisor = gcd_list(alkyl_protons)
    simplified_ratio = sorted([p // common_divisor for p in alkyl_protons])
    if simplified_ratio != sorted([2, 1, 1]):
        actual_ratio_str = ":".join(map(str, sorted(alkyl_protons)))
        return f"Incorrect. The alkyl signals for option {llm_answer} have an integration ratio of {actual_ratio_str}, which does not simplify to 2:1:1."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
print(check_correctness())