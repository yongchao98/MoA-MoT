import math
from functools import reduce

def get_gcd_for_list(numbers):
    """Helper function to calculate the greatest common divisor for a list of numbers."""
    if not numbers:
        return 1
    return reduce(math.gcd, numbers)

def check_correctness():
    """
    This function checks the correctness of the answer to the NMR spectroscopy question.
    It models the NMR properties of each compound and evaluates each mixture option
    against the given constraints.
    """
    # Define the predicted NMR spectra for each compound based on its structure and symmetry.
    # Each signal is a tuple of (integration, multiplicity).
    compounds = {
        "1,2,4,5-tetramethylbenzene": {
            "aromatic_signals": [(2, "singlet")],
            "alkyl_signals": [(12, "singlet")]
        },
        "1,2,3,5-tetramethylbenzene": {
            "aromatic_signals": [(1, "singlet"), (1, "singlet")],
            "alkyl_signals": [(6, "singlet"), (3, "singlet"), (3, "singlet")]
        },
        "1,2,3,4-tetramethylbenzene": {
            "aromatic_signals": [(2, "singlet")],
            "alkyl_signals": [(6, "singlet"), (6, "singlet")]
        },
        "1,4-diethylbenzene": {
            "aromatic_signals": [(4, "singlet")],
            "alkyl_signals": [(4, "quartet"), (6, "triplet")]
        }
    }

    # Define the options as presented in the question
    options = {
        "A": ("1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"),
        "B": ("1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"),
        "C": ("1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"),
        "D": ("1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene")
    }

    # Define the constraints from the question
    constraints = {
        "aromatic": {"count": 2, "multiplicity": "singlet", "ratio": (1, 1)},
        "alkyl": {"count": 3, "multiplicity": "singlet", "ratio": (2, 1, 1)}
    }
    
    provided_answer = "D"
    correct_option = None

    for option_key, (comp1_name, comp2_name) in options.items():
        comp1 = compounds[comp1_name]
        comp2 = compounds[comp2_name]

        # Combine the signals for the 1:1 mixture
        mixture_aromatic_signals = comp1["aromatic_signals"] + comp2["aromatic_signals"]
        mixture_alkyl_signals = comp1["alkyl_signals"] + comp2["alkyl_signals"]

        # --- Check Aromatic Constraints ---
        aromatic_ok = True
        # Check signal count
        if len(mixture_aromatic_signals) != constraints["aromatic"]["count"]:
            aromatic_ok = False
        # Check multiplicity
        if not all(sig[1] == constraints["aromatic"]["multiplicity"] for sig in mixture_aromatic_signals):
            aromatic_ok = False
        # Check ratio
        aromatic_integrations = sorted([sig[0] for sig in mixture_aromatic_signals])
        if aromatic_integrations[0] != aromatic_integrations[1]: # For 1:1 ratio, integrations must be equal
            aromatic_ok = False

        # --- Check Alkyl Constraints ---
        alkyl_ok = True
        # Check signal count
        if len(mixture_alkyl_signals) != constraints["alkyl"]["count"]:
            alkyl_ok = False
        # Check multiplicity
        if not all(sig[1] == constraints["alkyl"]["multiplicity"] for sig in mixture_alkyl_signals):
            alkyl_ok = False
        
        # Check ratio (only if other alkyl checks passed)
        if alkyl_ok:
            alkyl_integrations = sorted([sig[0] for sig in mixture_alkyl_signals], reverse=True)
            gcd = get_gcd_for_list(alkyl_integrations)
            simplified_ratio = tuple(i // gcd for i in alkyl_integrations)
            if simplified_ratio != constraints["alkyl"]["ratio"]:
                alkyl_ok = False

        if aromatic_ok and alkyl_ok:
            correct_option = option_key
            break # Found the correct option

    if correct_option == provided_answer:
        return "Correct"
    else:
        # Provide a detailed reason for the failure of the provided answer
        comp1_name, comp2_name = options[provided_answer]
        comp1 = compounds[comp1_name]
        comp2 = compounds[comp2_name]
        mixture_aromatic_signals = comp1["aromatic_signals"] + comp2["aromatic_signals"]
        mixture_alkyl_signals = comp1["alkyl_signals"] + comp2["alkyl_signals"]
        
        reason = f"The provided answer '{provided_answer}' is incorrect. The correct answer is '{correct_option}'.\n"
        reason += f"Let's re-evaluate the mixture for option {provided_answer} ({comp1_name} and {comp2_name}):\n"
        
        # Aromatic check for provided answer
        if len(mixture_aromatic_signals) != constraints["aromatic"]["count"]:
            reason += f"- Aromatic constraint failed: Expected {constraints['aromatic']['count']} signals, but the mixture produces {len(mixture_aromatic_signals)}.\n"
        
        # Alkyl check for provided answer
        if not all(sig[1] == constraints["alkyl"]["multiplicity"] for sig in mixture_alkyl_signals):
             reason += f"- Alkyl constraint failed: Expected all signals to be '{constraints['alkyl']['multiplicity']}', but the mixture contains other multiplicities.\n"
        elif len(mixture_alkyl_signals) != constraints["alkyl"]["count"]:
            reason += f"- Alkyl constraint failed: Expected {constraints['alkyl']['count']} signals, but the mixture produces {len(mixture_alkyl_signals)}.\n"
        else:
            alkyl_integrations = sorted([sig[0] for sig in mixture_alkyl_signals], reverse=True)
            gcd = get_gcd_for_list(alkyl_integrations)
            simplified_ratio = tuple(i // gcd for i in alkyl_integrations)
            if simplified_ratio != constraints["alkyl"]["ratio"]:
                 reason += f"- Alkyl constraint failed: Expected ratio {constraints['alkyl']['ratio']}, but the mixture produces a ratio of {simplified_ratio}.\n"

        return reason

# Run the check
print(check_correctness())