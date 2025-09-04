import re

def check_chemistry_answer():
    """
    Checks the correctness of the answer for the epoxide ring-opening reaction.
    
    This function verifies the two key constraints of the reaction:
    1. Regioselectivity: Attack at the less hindered carbon.
    2. Stereoselectivity: Inversion at the attacked center and retention at others.
    """
    
    # --- Problem Definition ---
    # Reactant: (1R,3R,4R,6S)-1,3,4-trimethyl-7-oxabicyclo[4.1.0]heptane
    # The epoxide is between C1 and C6.
    # C1 is quaternary (methyl-substituted), making it more hindered.
    # C6 is tertiary, making it less hindered.
    
    # The final answer from the LLM to be checked.
    llm_answer_option = "C"
    llm_answer_name = "(1R,2S,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol"

    # --- Verification Step 1: Regioselectivity ---
    # Attack must occur at the less hindered C6. This breaks the C6-O bond.
    # The -OH group forms at C1, and the new methyl group is added to C6.
    # When re-numbered by IUPAC rules, the C1(OH) is the new C1, and the adjacent
    # C6 (with the new methyl group) is the new C2.
    # Therefore, the product skeleton must be '1,2,4,5-tetramethylcyclohexan-1-ol'.
    # An incorrect attack at C1 would yield a '2,2,4,5-...' skeleton.
    
    expected_skeleton = "1,2,4,5-tetramethylcyclohexan-1-ol"
    if expected_skeleton not in llm_answer_name:
        return (f"Incorrect Regioselectivity: The product skeleton is wrong. "
                f"The nucleophile should attack the less hindered carbon (C6), leading to a '{expected_skeleton}' structure. "
                f"The answer's structure implies an incorrect attack at the more hindered C1.")

    # --- Verification Step 2: Stereochemistry ---
    
    # 2a. Parse the stereodescriptors from the answer's IUPAC name.
    match = re.search(r'\((.*?)\)', llm_answer_name)
    if not match:
        return "Stereochemistry Check Failed: Could not parse stereodescriptors from the answer name."
    
    descriptors_str = match.group(1)
    try:
        # Creates a dictionary like {1: 'R', 2: 'S', 4: 'R', 5: 'R'}
        descriptors = {int(s.strip()[:-1]): s.strip()[-1].upper() for s in descriptors_str.split(',')}
    except (ValueError, IndexError):
        return f"Stereochemistry Check Failed: Could not parse the descriptor string '{descriptors_str}'."

    # 2b. Check the stereocenters that are NOT attacked (retention of configuration).
    # Reactant centers: 1R, 3R, 4R
    # In the product, these correspond to C1, C5, and C4 respectively.
    # Expected retained configurations: 1R, 4R, 5R
    
    expected_retained = {1: 'R', 4: 'R', 5: 'R'}
    for center, config in expected_retained.items():
        if center not in descriptors or descriptors[center] != config:
            return (f"Stereochemistry Check Failed: The configuration at C{center} should be retained as '{config}' "
                    f"(from the original reactant), but the answer gives '{descriptors.get(center)}'.")

    # 2c. Check the stereocenter that IS attacked (inversion of configuration).
    # The attacked center is C6 in the reactant, which has an (S) configuration.
    # This becomes C2 in the product.
    # CRITICAL POINT: Geometric inversion (Sₙ2) does not guarantee the R/S label will flip.
    # The label depends on Cahn-Ingold-Prelog (CIP) priorities, which change during this reaction.
    # A rigorous 3D analysis shows that the geometric inversion at the C6(S) center
    # results in a product with a C2(S) configuration.
    
    expected_inverted_center_num = 2
    expected_inverted_center_config = 'S'
    
    actual_inverted_config = descriptors.get(expected_inverted_center_num)

    if actual_inverted_config != expected_inverted_center_config:
        # This part specifically checks for the most common error in this problem.
        if actual_inverted_config == 'R':
            return (f"Stereochemistry Check Failed: The configuration at the attacked center (C2 in the product) is incorrect. "
                    f"The answer assumes that inversion of a (S) center always results in a (R) center. This is a common error. "
                    f"Due to changes in substituent priorities (Cahn-Ingold-Prelog rules) during the reaction, the correct configuration for C2 is actually (S).")
        else:
            return (f"Stereochemistry Check Failed: The configuration at the attacked center (C2 in the product) should be '{expected_inverted_center_config}', "
                    f"but the answer gives '{actual_inverted_config}'.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_chemistry_answer()
print(result)