import re

def check_diels_alder_noesy_answer():
    """
    Checks the correctness of the given answer for the Diels-Alder NOESY problem.

    The function works by:
    1.  Defining the chemical principles of the reaction (Diels-Alder, endo-selectivity).
    2.  Identifying the specific protons that are spatially close in the major (endo) product.
    3.  Predicting the ¹H NMR signals (integration, multiplicity, chemical shift) for these interacting protons.
    4.  Parsing the signals described in the provided answer (Option D).
    5.  Comparing the parsed signals from the answer with the predicted signals of the interacting protons.
    6.  Returning "Correct" if they match, or a detailed reason if they do not.
    """

    # --- Step 1: Define Chemical Principles and Predictions ---

    # The reaction is a Diels-Alder between 1,2,3,4-tetramethyl-1,3-cyclopentadiene and maleic anhydride.
    # The major product is the 'endo' adduct due to kinetic control (secondary orbital overlap).
    
    # In the endo adduct, the key spatial proximity for a unique NOESY cross-peak is between:
    # Group 1: The two protons originating from the maleic anhydride.
    # Group 2: The six protons of the two methyl groups on the product's double bond.

    # Predict the ¹H NMR signals for these two interacting groups:
    # Group 1 (Anhydride Protons):
    # - Integration: 2H (two equivalent protons).
    # - Multiplicity: singlet (no adjacent protons to couple with).
    # - Chemical Shift: ~3.5 ppm (deshielded by adjacent carbonyls).
    predicted_interacting_signal_1 = {"integration": 2, "multiplicity": "singlet", "shift": 3.5}

    # Group 2 (Vinylic Methyl Protons):
    # - Integration: 6H (two equivalent methyl groups).
    # - Multiplicity: singlet (no coupling).
    # - Chemical Shift: ~1.7 ppm (typical for vinylic methyls).
    predicted_interacting_signal_2 = {"integration": 6, "multiplicity": "singlet", "shift": 1.7}

    # --- Step 2: Define and Parse the LLM's Answer ---

    llm_answer_choice = "D"
    options = {
        "A": "A 6H singlet at ~1 ppm and a 6H singlet at ~1.7 ppm",
        "B": "A 1H doublet at ~1.5 ppm and a 2H singlet at ~3.5 ppm",
        "C": "A 6H singlet at ~1 ppm and a 1H doublet at ~1.5 ppm",
        "D": "A 6H singlet at ~1.7 ppm and a 2H singlet at ~3.5 ppm"
    }
    llm_answer_text = options[llm_answer_choice]

    def parse_signal_description(desc):
        """Parses a string like 'A 6H singlet at ~1.7 ppm' into a dictionary."""
        # Regex to find integration, multiplicity, and shift
        match = re.search(r'(\d+)H\s+(singlet|doublet)\s+at\s+~?([\d\.]+)\s*ppm', desc)
        if not match:
            return None
        return {
            "integration": int(match.group(1)),
            "multiplicity": match.group(2),
            "shift": float(match.group(3))
        }

    # Split the answer text into two signal descriptions
    try:
        desc1, desc2 = llm_answer_text.split(' and ')
        parsed_signal_A = parse_signal_description(desc1)
        parsed_signal_B = parse_signal_description(desc2)
        if not parsed_signal_A or not parsed_signal_B:
            raise ValueError("Parsing failed")
    except (ValueError, AttributeError):
        return f"Error: Could not parse the signal descriptions in the answer: '{llm_answer_text}'"

    # --- Step 3: Check Correctness ---

    # Function to check if a parsed signal matches a predicted signal profile
    def signals_match(parsed, predicted):
        return (parsed['integration'] == predicted['integration'] and
                parsed['multiplicity'] == predicted['multiplicity'] and
                # Allow a small tolerance for chemical shift values
                abs(parsed['shift'] - predicted['shift']) < 0.2)

    # Check if the two parsed signals from the answer match the two predicted interacting signals.
    # The order of signals in the answer does not matter.
    is_correct_pair = (signals_match(parsed_signal_A, predicted_interacting_signal_1) and
                       signals_match(parsed_signal_B, predicted_interacting_signal_2)) or \
                      (signals_match(parsed_signal_A, predicted_interacting_signal_2) and
                       signals_match(parsed_signal_B, predicted_interacting_signal_1))

    if is_correct_pair:
        return "Correct"
    else:
        # If the answer is incorrect, provide a specific reason.
        # Check for signals that are inconsistent with the product structure.
        if parsed_signal_A['multiplicity'] == 'doublet' or parsed_signal_B['multiplicity'] == 'doublet':
            return ("Incorrect. The answer includes a 'doublet'. The major product is highly symmetrical, "
                    "and the key proton groups (methyls, anhydride protons) are expected to be singlets. "
                    "A '1H doublet' is inconsistent with the product's structure.")

        # Check if the answer pairs two signals that exist but are not spatially close (e.g., Option A).
        bridgehead_methyl_signal = {"integration": 6, "multiplicity": "singlet", "shift": 1.0}
        is_bridgehead_A = signals_match(parsed_signal_A, bridgehead_methyl_signal)
        is_bridgehead_B = signals_match(parsed_signal_B, bridgehead_methyl_signal)
        is_vinylic_A = signals_match(parsed_signal_A, predicted_interacting_signal_2)
        is_vinylic_B = signals_match(parsed_signal_B, predicted_interacting_signal_2)

        if (is_bridgehead_A and is_vinylic_B) or (is_bridgehead_B and is_vinylic_A):
            return ("Incorrect. The answer pairs the vinylic methyls (~1.7 ppm) with the bridgehead methyls (~1.0 ppm). "
                    "While both signals are expected in the spectrum, these two groups are spatially distant in the major (endo) product "
                    "and would not show a strong NOESY cross-peak. The cross-peak arises from the proximity of the vinylic methyls "
                    "to the anhydride protons (~3.5 ppm).")

        # Generic failure message.
        return ("Incorrect. The answer does not correctly identify the pair of signals corresponding to the protons that are "
                "spatially close in the major product. The NOESY cross-peak is expected between the vinylic methyl protons "
                "(6H singlet at ~1.7 ppm) and the anhydride protons (2H singlet at ~3.5 ppm).")

# Execute the checker function and print the result.
result = check_diels_alder_noesy_answer()
# print(result) # This would print "Correct"