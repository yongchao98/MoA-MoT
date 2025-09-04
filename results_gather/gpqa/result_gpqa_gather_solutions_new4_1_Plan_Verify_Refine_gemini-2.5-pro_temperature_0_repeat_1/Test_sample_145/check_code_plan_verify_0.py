import re

def check_diels_alder_product():
    """
    Checks the correctness of the answer for the Diels-Alder reaction between
    5-fluorocyclopenta-1,3-diene and maleic anhydride.

    The check is based on established stereochemical principles:
    1. The Alder-Endo Rule for kinetic control.
    2. The principle of syn-facial selectivity for C5-electronegative substituents.
    3. The resulting product geometry and its corresponding IUPAC nomenclature.
    """
    
    # The options provided in the problem
    options = {
        'A': "(3aR,4R,7S,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        'B': "(3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        'C': "(3aR,4S,7R,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        'D': "(3aR,4R,7S,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione"
    }
    
    # The final answer to be checked
    proposed_answer_key = 'C'
    
    # --- Define Chemical Principles and Nomenclature Mappings ---

    # Principle 1: The major product is the 'endo' adduct.
    # Mapping: The 'endo' adduct has a (4S, 7R) or (4R, 7S) configuration depending on the enantiomer.
    # We will check for the relative stereochemistry that defines the endo adduct.
    # Based on literature, the (...,4S,7R,...) and (...,4R,7S,...) descriptors define the endo/exo skeletons.
    # The consensus is that (...,4S,7R,...) corresponds to the endo adduct for the (3aR, 7aS) enantiomer.
    endo_skeleton_descriptors = ('4S', '7R')
    
    # Principle 2: Facial selectivity for C5-F is 'syn-attack' due to electronic effects.
    # This 'syn-attack' leads to a product where the F is 'anti' to the anhydride ring.
    # This 'anti-F' geometry in an 'endo' framework corresponds to the '8s' descriptor.
    predicted_c8_descriptor = '8s'

    # --- Analysis of the Proposed Answer ---

    chosen_answer_name = options.get(proposed_answer_key)
    if not chosen_answer_name:
        return f"Invalid answer key '{proposed_answer_key}' provided."

    # Use regex to parse the stereodescriptors from the IUPAC name
    match = re.search(r'\((.*?)\)', chosen_answer_name)
    if not match:
        return f"Could not parse stereodescriptors from the name of option {proposed_answer_key}."
    
    descriptors_str = match.group(1)
    descriptors = [d.strip() for d in descriptors_str.split(',')]
    
    try:
        # Extract the descriptors needed for the check
        actual_c4_desc = descriptors[1]
        actual_c7_desc = descriptors[2]
        actual_c8_desc = descriptors[4]
    except IndexError:
        return f"Could not extract all necessary stereodescriptors from option {proposed_answer_key}."

    # --- Verification against Principles ---

    # Check 1: Endo/Exo Selectivity
    # The major product should be the 'endo' adduct.
    actual_skeleton_descriptors = (actual_c4_desc, actual_c7_desc)
    if actual_skeleton_descriptors != endo_skeleton_descriptors:
        return (f"Incorrect. The answer {proposed_answer_key} violates the Endo Rule. "
                f"The major product should be the 'endo' adduct, which has the stereochemical signature "
                f"{endo_skeleton_descriptors} for carbons (C4, C7). "
                f"Option {proposed_answer_key} has {actual_skeleton_descriptors}, which corresponds to the 'exo' adduct.")

    # Check 2: Facial Selectivity and C8 Stereochemistry
    # The major product results from syn-facial attack, leading to an anti-F geometry, which is '8s'.
    if actual_c8_desc != predicted_c8_descriptor:
        return (f"Incorrect. The answer {proposed_answer_key} has the wrong stereochemistry for the fluorine substituent. "
                f"The major pathway is syn-facial attack, which leads to a product with an 'anti' fluorine. "
                f"This geometry corresponds to the C8 descriptor '{predicted_c8_descriptor}'. "
                f"However, option {proposed_answer_key} has the descriptor '{actual_c8_desc}'.")

    # If all checks pass, the answer is consistent with the established principles.
    return "Correct"

# Execute the check and print the result
result = check_diels_alder_product()
print(result)