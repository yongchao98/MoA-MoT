import re

def check_correctness_of_answer(llm_choice):
    """
    This function checks the correctness of a given answer choice (A, B, C, or D)
    based on the constraints from the chemistry problem.

    The problem describes a di-substituted 6-membered aromatic ring with 8 total carbons,
    a carbonyl group, and an aromatic-halogen bond.

    This implies the structure is an acetyl-halo-benzene (e.g., 4-chloroacetophenone),
    with the molecular formula C8H7XO.
    """

    # --- Step 1: Define constraints based on the deduced structure C8H7XO ---
    expected_total_protons = 7
    expected_aromatic_protons = 4  # A di-substituted benzene ring has 4 protons.
    
    # --- Step 2: Define the NMR data for each option ---
    options = {
        "A": "1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)",
        "B": "1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)",
        "C": "1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)",
        "D": "1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)"
    }

    if llm_choice not in options:
        return f"Invalid choice '{llm_choice}'. Please select from A, B, C, or D."

    nmr_string = options[llm_choice]

    # --- Step 3: Parse the NMR data string for the chosen option ---
    try:
        # Regex to find patterns like "7.8 (2H, d)"
        signals = re.findall(r"(\d+\.?\d*)\s*\((\d+)H,\s*([a-z]+)\)", nmr_string)
        parsed_signals = [{
            "shift": float(shift),
            "protons": int(protons),
            "splitting": splitting
        } for shift, protons, splitting in signals]
    except (ValueError, IndexError):
        return f"Error: Could not parse the NMR data for option {llm_choice}."

    # --- Step 4: Check the parsed data against the constraints ---

    # Constraint Check 1: Total number of protons
    total_protons = sum(s['protons'] for s in parsed_signals)
    if total_protons != expected_total_protons:
        return (f"Incorrect. The total number of protons is {total_protons}, but the expected "
                f"structure C8H7XO should have {expected_total_protons} protons.")

    # Constraint Check 2: Number of aromatic protons
    aromatic_protons = sum(s['protons'] for s in parsed_signals if 6.5 <= s['shift'] <= 8.5)
    if aromatic_protons != expected_aromatic_protons:
        return (f"Incorrect. The number of aromatic protons is {aromatic_protons}, but a "
                f"di-substituted benzene ring must have {expected_aromatic_protons}.")

    # Constraint Check 3: Presence of an acetyl group signal (-COCH3)
    # Expected: 3 protons, singlet splitting, chemical shift ~2.0-2.8 ppm
    has_acetyl_signal = any(
        s['protons'] == 3 and s['splitting'] == 's' and 2.0 <= s['shift'] <= 2.8
        for s in parsed_signals
    )
    if not has_acetyl_signal:
        return ("Incorrect. The NMR data does not contain a signal consistent with an acetyl group "
                "(which should be a 3H singlet around 2.3 ppm).")

    # If all constraints are met, the answer is correct.
    return "Correct"

# Based on the analysis, the correct option is 'A'.
# We will run the check on 'A' to confirm.
result = check_correctness_of_answer("A")
print(result)

# To be thorough, let's see why the other options are wrong.
# print("\nChecking other options for failure reasons:")
# print(f"B: {check_correctness_of_answer('B')}")
# print(f"C: {check_correctness_of_answer('C')}")
# print(f"D: {check_correctness_of_answer('D')}")
