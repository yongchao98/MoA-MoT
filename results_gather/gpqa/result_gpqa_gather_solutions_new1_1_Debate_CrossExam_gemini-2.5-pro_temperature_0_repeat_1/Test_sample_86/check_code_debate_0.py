import re

def check_nmr_answer(final_answer_str):
    """
    Checks the correctness of the selected 1H NMR data based on the problem description.

    The problem describes a di-substituted 6-membered aromatic ring with 8 total carbons,
    a carbonyl group, and an aromatic-halogen bond. This uniquely points to a haloacetophenone
    structure (X-C6H4-COCH3). This script checks which NMR data option matches this structure.
    """

    # --- Step 1: Define the problem's options and the expected properties ---

    options = {
        "A": "1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)",
        "B": "1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)",
        "C": "1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)",
        "D": "1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)"
    }

    # Properties deduced from the structure (haloacetophenone, X-C6H4-COCH3):
    expected_total_protons = 7  # 4 on the aromatic ring, 3 on the methyl group.
    expected_aromatic_protons = 4

    # --- Step 2: Create a function to parse and validate a spectrum ---

    def parse_and_validate(spectrum_str):
        """Parses an NMR data string and checks it against chemical rules."""
        
        # Use regex to find all signals in the format: shift (integrationH, multiplicity)
        pattern = re.compile(r"(\d+\.?\d*)\s*\((\d+)H,\s*([a-z]+)\)")
        signals = pattern.findall(spectrum_str)
        
        if not signals:
            return f"Could not parse any valid NMR signals from the string."

        parsed_signals = [{
            "shift": float(shift),
            "integration": int(integration),
            "multiplicity": multiplicity
        } for shift, integration, multiplicity in signals]

        # Constraint 1: Check total number of protons.
        total_protons_found = sum(s['integration'] for s in parsed_signals)
        if total_protons_found != expected_total_protons:
            return f"Incorrect total proton count. Expected {expected_total_protons}, but found {total_protons_found}."

        # Constraint 2: Check for aromatic protons.
        aromatic_signals = [s for s in parsed_signals if 6.5 <= s['shift'] <= 8.5]
        aromatic_protons_found = sum(s['integration'] for s in aromatic_signals)
        if aromatic_protons_found != expected_aromatic_protons:
            return f"Incorrect number of aromatic protons. Expected {expected_aromatic_protons}, but found {aromatic_protons_found} in the 6.5-8.5 ppm range."

        # Constraint 3: Check for the acetyl methyl group signal.
        # Must be a 3H singlet in the typical methyl ketone range.
        acetyl_signal = [s for s in parsed_signals if 2.0 <= s['shift'] <= 2.7 and s['integration'] == 3 and s['multiplicity'] == 's']
        if not acetyl_signal:
            return "Missing the characteristic acetyl methyl signal (a 3H singlet between 2.0-2.7 ppm)."

        # Constraint 4: Check for an aldehyde proton (which would indicate an incorrect structure).
        aldehyde_signal = [s for s in parsed_signals if s['shift'] >= 9.0]
        if aldehyde_signal:
            return "An aldehyde proton signal was found (~9.0-10.0 ppm), which is inconsistent with the deduced ketone structure."
            
        # Constraint 5: Check for the para-substitution pattern (a strong positive indicator).
        # The aromatic region should contain two 2H doublets.
        para_pattern_signals = [s for s in aromatic_signals if s['integration'] == 2 and s['multiplicity'] == 'd']
        if len(para_pattern_signals) != 2:
             return "The aromatic signals do not match the classic para-substitution pattern (two 2H doublets)."

        # If all checks pass, the spectrum is valid.
        return None

    # --- Step 3: Execute the check on the provided answer ---

    # Extract the answer key (e.g., 'B' from '<<<B>>>')
    try:
        answer_key = re.search(r'<<<([A-D])>>>', final_answer_str).group(1)
    except (AttributeError, IndexError):
        return f"Invalid answer format: {final_answer_str}. Expected format like '<<<A>>>'."

    if answer_key not in options:
        return f"Invalid answer key '{answer_key}'. Must be one of {list(options.keys())}."

    # Validate the spectrum corresponding to the answer key
    validation_result = parse_and_validate(options[answer_key])

    if validation_result is None:
        return "Correct"
    else:
        return f"Incorrect. The data for option {answer_key} is not consistent with the deduced structure (para-haloacetophenone). Reason: {validation_result}"

# The final answer provided in the prompt to be checked
final_answer_from_llm = "<<<B>>>"

# Run the checking code and print the result
result = check_nmr_answer(final_answer_from_llm)
print(result)