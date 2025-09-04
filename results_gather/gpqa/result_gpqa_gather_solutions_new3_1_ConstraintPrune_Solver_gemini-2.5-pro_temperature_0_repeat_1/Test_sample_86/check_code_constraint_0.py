import re

def parse_nmr_data(data_string):
    """Parses a string of NMR data into a list of signal dictionaries."""
    signals = []
    # Regex to find patterns like: 7.8 (2H, d)
    # It captures the chemical shift, integration, and splitting
    pattern = re.compile(r"(\d+\.?\d*)\s*\((\d+)H,\s*([a-z]+)\)")
    matches = pattern.findall(data_string)
    for match in matches:
        signals.append({
            "shift": float(match[0]),
            "integration": int(match[1]),
            "splitting": match[2]
        })
    return signals

def check_answer(question_details, options, final_answer_letter):
    """
    Checks the correctness of the final answer based on the question's constraints.

    Args:
        question_details (dict): A dictionary containing the constraints from the question.
        options (dict): A dictionary of the NMR data for each option letter.
        final_answer_letter (str): The letter of the chosen final answer (e.g., 'C').

    Returns:
        str: "Correct" if the answer is correct, otherwise a reason for the failure.
    """
    # Step 1: Deduce the expected structure from the question's constraints.
    # Aromatic ring = 6 carbons. Total carbons = 8. -> Substituents have 2 carbons.
    # Substituents = Halogen (0 carbons) + a group with 2 carbons and a C=O.
    # This second group must be an acetyl group (-COCH3).
    # So, the structure is a haloacetophenone.
    # The total number of hydrogens should be 7 (from C8H7XO).
    
    expected_structure = "para-haloacetophenone"
    
    # Step 2: Predict the 1H NMR spectrum for the expected structure (para-isomer is most likely).
    # - Aromatic region: Two doublets, each for 2H, in the 7.0-8.5 ppm range.
    # - Acetyl group region: One singlet for 3H, in the 2.1-2.7 ppm range.
    
    # Step 3: Analyze the chosen answer (Option C).
    if final_answer_letter not in options:
        return f"Invalid answer letter '{final_answer_letter}'. Valid options are {list(options.keys())}."
        
    chosen_option_data = options[final_answer_letter]
    signals = parse_nmr_data(chosen_option_data)
    
    # Check total proton count
    total_protons = sum(s['integration'] for s in signals)
    if total_protons != 7:
        return (f"Incorrect. The total proton integration for option {final_answer_letter} is {total_protons}, "
                f"but the expected structure (C8H7XO) should have 7 protons.")

    # Check for aromatic signals
    aromatic_signals = [s for s in signals if 6.5 <= s['shift'] <= 8.5]
    if len(aromatic_signals) != 2:
        return (f"Incorrect. Option {final_answer_letter} has {len(aromatic_signals)} signals in the aromatic region. "
                f"Expected 2 signals for a para-substituted ring.")
    
    aromatic_doublets_2H = [s for s in aromatic_signals if s['splitting'] == 'd' and s['integration'] == 2]
    if len(aromatic_doublets_2H) != 2:
        return (f"Incorrect. The aromatic signals in option {final_answer_letter} do not match the expected pattern "
                f"of two 2H doublets for a para-substituted ring.")

    # Check for acetyl group signal
    acetyl_signals = [s for s in signals if 2.0 <= s['shift'] <= 2.7]
    if len(acetyl_signals) != 1:
        return (f"Incorrect. Option {final_answer_letter} has {len(acetyl_signals)} signals in the typical acetyl region. "
                f"Expected exactly one.")
        
    acetyl_singlet_3H = [s for s in acetyl_signals if s['splitting'] == 's' and s['integration'] == 3]
    if len(acetyl_singlet_3H) != 1:
        return (f"Incorrect. The signal in the acetyl region in option {final_answer_letter} is not a 3H singlet as expected.")

    # Check for extraneous signals (like an aldehyde proton)
    aldehyde_signals = [s for s in signals if s['shift'] > 9.0]
    if aldehyde_signals:
        return (f"Incorrect. Option {final_answer_letter} contains a signal at {aldehyde_signals[0]['shift']} ppm, "
                f"which is characteristic of an aldehyde, not the expected ketone.")

    # If all checks for the chosen option pass, it is correct.
    # We can also verify why other options are wrong.
    for letter, data in options.items():
        if letter == final_answer_letter:
            continue
        # Simple checks for other options
        other_signals = parse_nmr_data(data)
        if not any(6.5 <= s['shift'] <= 8.5 for s in other_signals) and letter == 'D':
             # Option D is wrong because it has no aromatic signals
             pass
        if any(s['shift'] > 9.0 for s in other_signals) and letter == 'A':
             # Option A is wrong because it has an aldehyde signal
             pass
        if len([s for s in other_signals if 6.5 <= s['shift'] <= 8.5]) != 4 and letter == 'B':
             # Option B is wrong because it doesn't have 4 aromatic protons
             pass

    return "Correct"

# --- Main execution ---
# Define the problem constraints and options from the prompt
question_constraints = {
    "carbons": 8,
    "structure": "di-substituted 6-membered aromatic ring",
    "functional_groups": ["carbonyl", "aromatic-halogen"]
}

# Options as presented in the original question
options = {
    "A": "1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)",
    "B": "1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)",
    "C": "1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)",
    "D": "1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)"
}

# The final answer to be checked
final_answer = "C"

# Run the check
result = check_answer(question_constraints, options, final_answer)
print(result)