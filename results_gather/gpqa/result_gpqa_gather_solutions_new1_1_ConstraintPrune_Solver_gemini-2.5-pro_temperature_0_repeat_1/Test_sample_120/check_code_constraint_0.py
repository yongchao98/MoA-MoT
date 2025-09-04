import re

def check_answer():
    """
    Checks the correctness of the answer to the organic chemistry question.
    """
    # --- Problem Definition & Analysis ---

    # 1. Regioselectivity Analysis:
    # The nucleophile (Me-) attacks the less hindered epoxide carbon.
    # C1 is quaternary (more hindered), C6 is tertiary (less hindered).
    # Attack must occur at C6.
    # This results in a product with a 1,2,4,5-tetramethylcyclohexan-1-ol skeleton.
    expected_skeleton = "1,2,4,5-tetramethylcyclohexan-1-ol"

    # 2. Stereochemistry Analysis:
    # The reaction is an SN2 attack, causing inversion at the attacked center.
    # Original configurations: (1R, 3R, 4R, 6S)
    # - Attack at C6 (S) causes inversion -> new C2 becomes (R).
    # - C1 (R) is retained -> new C1 remains (R).
    # - C3 (R) is retained -> new C5 remains (R).
    # - C4 (R) is retained -> new C4 remains (R).
    # The expected stereochemistry is (1R, 2R, 4R, 5R).
    expected_stereochem_str = "(1R,2R,4R,5R)"

    # The expected full name of the product.
    expected_product_name = f"{expected_stereochem_str}-{expected_skeleton}"

    # --- The Answer to Check ---
    # The final answer from the LLM is C.
    # Let's get the text for option C from the question prompt.
    option_c_text = "(1R,2R,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol"
    
    # --- Verification ---

    # Check 1: Regioselectivity (Product Skeleton)
    # The skeleton is the part of the name describing the carbon framework and functional groups.
    if expected_skeleton not in option_c_text:
        # This would happen if the answer chose a product from attack at C1.
        return (f"Incorrect. The regioselectivity is wrong. The reaction should yield a "
                f"'{expected_skeleton}' skeleton by attacking the less hindered C6, but the "
                f"answer implies a different skeleton.")

    # Check 2: Stereochemistry
    # Extract the stereochemical descriptors from the answer.
    match = re.match(r'(\(.*\))', option_c_text)
    if not match:
        return f"Incorrect. Could not parse stereochemistry from the answer '{option_c_text}'."
    
    actual_stereochem_str = match.group(1)

    if actual_stereochem_str != expected_stereochem_str:
        reasons = []
        if "1R" not in actual_stereochem_str:
            reasons.append("C1 should be R (retention)")
        if "2R" not in actual_stereochem_str:
            reasons.append("C2 should be R (inversion from S)")
        if "4R" not in actual_stereochem_str:
            reasons.append("C4 should be R (retention)")
        if "5R" not in actual_stereochem_str:
            reasons.append("C5 should be R (retention)")
        
        return (f"Incorrect. The stereochemistry is wrong. Expected {expected_stereochem_str} "
                f"but got {actual_stereochem_str}. Specific errors: {'; '.join(reasons)}.")

    # If both checks pass, the answer is correct.
    if option_c_text == expected_product_name:
        return "Correct"
    else:
        # Fallback for any other discrepancy.
        return f"Incorrect. The expected product is '{expected_product_name}', but the answer is '{option_c_text}'."

# Execute the check and print the result.
result = check_answer()
print(result)