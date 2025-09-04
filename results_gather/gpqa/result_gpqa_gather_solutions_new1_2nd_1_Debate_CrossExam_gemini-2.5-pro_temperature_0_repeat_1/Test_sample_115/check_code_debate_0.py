import re

def check_chemistry_answer():
    """
    This function checks the correctness of the provided answer for the multi-step synthesis problem.
    It verifies the two main constraints of the final Diels-Alder reaction:
    1. Product Skeleton (Regiochemistry)
    2. Product Stereochemistry (from the cis-dienophile)
    """

    # The four options provided in the original question
    options = {
        "A": "(1S,4S)-4,6,6-trimethylcyclohex-2-enol",
        "B": "(1S,4R,5S,6S)-4,5,6-trimethylcyclohex-2-enol",
        "C": "(1S,4R,5S,6R)-4,5,6-trimethylcyclohex-2-enol",
        "D": "(1S,4R)-4,6,6-trimethylcyclohex-2-enol"
    }

    # The final answer provided by the LLM being checked
    llm_answer = "C"

    # --- Constraint 1: Product Skeleton (Regiochemistry) ---
    # The Diels-Alder reaction is between cis-but-2-ene (a C4 dienophile) and
    # (1E,3E)-penta-1,3-dien-1-ol (a C5 diene).
    # The resulting product is a cyclohexene ring with three methyl groups and one alcohol.
    # IUPAC numbering places the substituents at positions 4, 5, and 6.
    # Therefore, the correct skeleton must be "4,5,6-trimethylcyclohex-2-enol".
    
    correct_skeleton_pattern = "4,5,6-trimethylcyclohex-2-enol"
    
    candidates_after_skeleton_check = []
    for option_key, option_name in options.items():
        if correct_skeleton_pattern in option_name:
            candidates_after_skeleton_check.append(option_key)

    # Check if the LLM's answer has the correct skeleton.
    if llm_answer not in candidates_after_skeleton_check:
        incorrect_skeleton = "4,6,6-trimethylcyclohex-2-enol" if "4,6,6" in options[llm_answer] else "an unknown skeleton"
        return (f"Incorrect. The provided answer '{llm_answer}' has the wrong product skeleton. "
                f"The Diels-Alder reaction should produce a '{correct_skeleton_pattern}' skeleton, "
                f"but option {llm_answer} has a '{incorrect_skeleton}' skeleton.")

    # --- Constraint 2: Product Stereochemistry ---
    # The Diels-Alder reaction is stereospecific. The dienophile is cis-but-2-ene.
    # This means the two methyl groups it provides (at C5 and C6) must be cis to each other in the product.
    # Rule for adjacent stereocenters:
    # - cis relationship corresponds to opposite R/S descriptors (R,S or S,R).
    # - trans relationship corresponds to identical R/S descriptors (R,R or S,S).
    
    final_candidate = None
    
    for option_key in candidates_after_skeleton_check:
        option_name = options[option_key]
        
        # Use regex to find the stereochemistry at C5 and C6
        match = re.search(r'5(S|R),6(S|R)', option_name)
        
        if not match:
            # This option has the right skeleton but lacks the necessary stereocenters to check.
            continue
            
        c5_stereo = match.group(1)
        c6_stereo = match.group(2)
        
        # Check for a cis relationship (opposite R/S descriptors)
        if c5_stereo != c6_stereo:
            # This option has the correct cis stereochemistry.
            if final_candidate is not None:
                # This would mean more than one option is correct, which is an issue with the question itself.
                return "Error in checking logic: Found multiple options satisfying all constraints."
            final_candidate = option_key

    # --- Final Verdict ---
    if final_candidate is None:
        return "Incorrect. No option satisfies all constraints. The provided answer '{llm_answer}' is incorrect because no option has the required 'cis' stereochemistry at C5 and C6."

    if final_candidate == llm_answer:
        return "Correct"
    else:
        # This means the LLM chose an answer that passed the skeleton check but failed the stereochem check.
        llm_option_name = options[llm_answer]
        match = re.search(r'5(S|R),6(S|R)', llm_option_name)
        c5 = match.group(1)
        c6 = match.group(2)
        relationship = "trans" # Since it failed the check, it must be trans
        
        return (f"Incorrect. The provided answer '{llm_answer}' is wrong. "
                f"While it has the correct skeleton, its stereochemistry at C5 and C6 is ({c5}, {c6}), "
                f"which corresponds to a '{relationship}' relationship. The reaction requires a 'cis' relationship "
                f"because the dienophile is cis-but-2-ene. The correct answer is '{final_candidate}'.")

# Execute the check and print the result
result = check_chemistry_answer()
print(result)