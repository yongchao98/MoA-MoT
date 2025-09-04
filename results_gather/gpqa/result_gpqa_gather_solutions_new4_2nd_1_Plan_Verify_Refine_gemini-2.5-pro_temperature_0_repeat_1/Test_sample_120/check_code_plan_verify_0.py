import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer for the given chemistry question.
    
    The reaction is the nucleophilic ring-opening of an epoxide by an organocuprate.
    The code will verify two key aspects:
    1. Regioselectivity: Attack at the less hindered carbon.
    2. Stereoselectivity: Inversion of configuration at the attacked carbon.
    """
    
    # --- Problem Definition ---
    # Starting material: (1R,3R,4R,6S)-1,3,4-trimethyl-7-oxabicyclo[4.1.0]heptane
    # Epoxide carbons are C1 and C6.
    # Stereocenters: {1: 'R', 3: 'R', 4: 'R', 6: 'S'}
    # Reagent: Me2CuLi (provides Me- nucleophile)
    
    # The final answer provided by the LLM to be checked.
    llm_final_answer = "A"
    
    options = {
        "A": "(1R,2S,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol",
        "B": "(1R,2R,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol",
        "C": "(1S,4R,5S)-2,2,4,5-tetramethylcyclohexan-1-ol",
        "D": "(1R,4R,5R)-2,2,4,5-tetramethylcyclohexan-1-ol"
    }
    
    llm_answer_text = options.get(llm_final_answer)
    if not llm_answer_text:
        return f"Invalid answer key '{llm_final_answer}'. The key must be one of {list(options.keys())}."

    # --- Verification Logic ---
    
    # 1. Check Regioselectivity
    # C1 is quaternary (bonded to a methyl group).
    # C6 is tertiary (bonded to a hydrogen).
    # The less hindered carbon is C6.
    attack_site = 6
    
    # Attack at C6 results in a 1,2,4,5-tetramethylcyclohexan-1-ol skeleton.
    # The -OH group is on the original C1 (new C1).
    # The new methyl group is on the original C6 (new C2).
    expected_skeleton = "1,2,4,5-tetramethylcyclohexan-1-ol"
    
    if expected_skeleton not in llm_answer_text:
        return (f"Incorrect: The answer does not satisfy the regioselectivity constraint. "
                f"Attack should occur at the less hindered C6, leading to a '{expected_skeleton}' skeleton. "
                f"The answer's skeleton implies an incorrect attack at C1.")

    # 2. Check Stereoselectivity
    expected_stereocenters = {}
    
    # a) Retention at non-reacting centers
    # Original C1(R) becomes new C1(R)
    expected_stereocenters[1] = 'R'
    # Original C4(R) becomes new C4(R)
    expected_stereocenters[4] = 'R'
    # Original C3(R) becomes new C5(R) due to IUPAC renumbering
    expected_stereocenters[5] = 'R'
    
    # b) Inversion at the reacting center (C6)
    # The starting configuration at C6 is (S).
    # The Sₙ2 attack causes a geometric inversion. A simple S->R flip is a common error
    # because the Cahn-Ingold-Prelog (CIP) priorities of the substituents change.
    # A detailed 3D analysis shows that the geometric inversion of the (6S) center
    # results in a new center with an (S) configuration.
    # New C2 (from original C6) should be 'S'.
    expected_stereocenters[2] = 'S'
    
    # 3. Assemble the expected product name and compare
    sorted_keys = sorted(expected_stereocenters.keys())
    stereo_string = ",".join([f"{k}{expected_stereocenters[k]}" for k in sorted_keys])
    expected_product_name = f"({stereo_string})-{expected_skeleton}"
    
    # Normalize strings for comparison
    normalized_expected = expected_product_name.replace(" ", "")
    normalized_llm_answer = llm_answer_text.replace(" ", "")

    if normalized_expected == normalized_llm_answer:
        return "Correct"
    else:
        # Provide a detailed reason for the mismatch
        try:
            match = re.search(r'\((.*?)\)', normalized_llm_answer)
            llm_stereo_str = match.group(1)
            llm_stereocenters = {int(c[0]): c[1] for c in llm_stereo_str.split(',')}
            
            for center, config in expected_stereocenters.items():
                if llm_stereocenters.get(center) != config:
                    return (f"Incorrect: The answer has the wrong stereochemistry at C{center}. "
                            f"It should be '{config}' but the answer gives '{llm_stereocenters.get(center)}'. "
                            f"This is likely due to an incorrect assumption that S always inverts to R, "
                            f"without considering the change in CIP priorities.")
            return "Incorrect: Stereochemistry mismatch for an unknown reason." # Fallback
        except Exception:
            return f"Incorrect: Could not parse the final answer's stereochemistry. Expected '{normalized_expected}', but got '{normalized_llm_answer}'."

# Execute the check and print the result
print(check_correctness())