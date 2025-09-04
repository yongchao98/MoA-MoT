import re

def check_organic_reaction_answer():
    """
    Checks the correctness of the answer for the given organocuprate reaction.
    This function codifies the chemical rules to derive the expected product and
    compares it with the provided answer.
    """
    # --- Problem Definition & Rules ---
    # Reactant: (1R,3R,4R,6S)-1,3,4-trimethyl-7-oxabicyclo[4.1.0]heptane
    # Rule 1: Attack at the less hindered epoxide carbon (C1 vs C6).
    # Rule 2: Inversion of configuration at the attacked carbon.
    
    # --- Step 1: Analyze Reactant and Site of Attack ---
    reactant_config = {'1': 'R', '3': 'R', '4': 'R', '6': 'S'}
    # C1 has a methyl group, C6 does not. C6 is less hindered.
    attacked_carbon_old_idx = '6'

    # --- Step 2: Determine Product Connectivity and Numbering ---
    # Attack at C6 adds a methyl group. Epoxide opens to form -OH at C1.
    # New IUPAC numbering: C-OH is C1. Number towards next substituent.
    # Old C1 -> New C1
    # Old C6 -> New C2
    # Old C4 -> New C4
    # Old C3 -> New C5
    expected_base_name = "1,2,4,5-tetramethylcyclohexan-1-ol"

    # --- Step 3: Predict Product Stereochemistry ---
    expected_product_config = {}

    # At New C1 (from Old C1): Configuration is retained as it's not the SN2 center.
    expected_product_config['1'] = reactant_config['1']  # Expected: R

    # At New C2 (from Old C6): Attacked carbon, configuration inverts.
    expected_product_config['2'] = 'R' if reactant_config[attacked_carbon_old_idx] == 'S' else 'S' # Expected: R

    # At New C4 (from Old C4): Spectator center, configuration retained.
    expected_product_config['4'] = reactant_config['4']  # Expected: R

    # At New C5 (from Old C3): Spectator center, configuration retained.
    expected_product_config['5'] = reactant_config['3']  # Expected: R

    # --- Step 4: Compare with the Provided Answer ---
    # The provided answer is C, which corresponds to:
    llm_answer_text = "(1R,2R,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol"

    # Parse the answer string
    match = re.match(r'\((.*?)\)-(.*)', llm_answer_text)
    if not match:
        return f"Could not parse the answer format: {llm_answer_text}"
    
    answer_stereo_part = match.group(1)
    answer_base_name = match.group(2)

    # Check 1: Base name (connectivity)
    if answer_base_name != expected_base_name:
        return (f"Incorrect molecular structure. The reaction should yield "
                f"'{expected_base_name}', but the answer gives '{answer_base_name}'.")

    # Check 2: Stereochemistry
    try:
        answer_stereo_pairs = answer_stereo_part.split(',')
        answer_config = {p[0]: p[1].upper() for p in answer_stereo_pairs}
    except Exception:
        return f"Could not parse the stereochemistry from the answer: {answer_stereo_part}"

    mismatches = []
    for center in sorted(expected_product_config.keys()):
        if expected_product_config.get(center) != answer_config.get(center):
            mismatches.append(
                f"At C{center}, expected configuration is {expected_product_config.get(center)} "
                f"but the answer gives {answer_config.get(center)}."
            )
            
    if mismatches:
        error_report = "Incorrect stereochemistry:\n" + "\n".join(mismatches)
        # Add reasoning for the first mismatch found
        first_mismatch_center = mismatches[0][5]
        if first_mismatch_center == '2':
            error_report += "\nReasoning: The nucleophile attacks the less hindered C6, causing inversion of its (S) configuration to (R). This carbon is C2 in the product, so it must be (2R)."
        elif first_mismatch_center == '1':
            error_report += "\nReasoning: The configuration at C1 is retained because it is not the site of nucleophilic attack. It was (1R) and should remain (1R)."
        else: # C4 or C5
            error_report += f"\nReasoning: The stereocenter at C{first_mismatch_center} is a spectator and its configuration should be retained from the reactant."
        return error_report

    return "Correct"

# Execute the check
result = check_organic_reaction_answer()
print(result)