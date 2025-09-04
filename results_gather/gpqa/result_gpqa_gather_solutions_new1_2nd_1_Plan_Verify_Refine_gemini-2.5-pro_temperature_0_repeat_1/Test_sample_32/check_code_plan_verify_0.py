import re

def check_answer_correctness():
    """
    This function checks the correctness of the provided answer for the given chemistry question.
    It verifies the answer against two main constraints derived from the problem statement:
    1. Structural Correctness: The product name must match the reactants (thiophene -> epithio bridge).
    2. Stereochemical Correctness: The product must be the EXO isomer, which has a specific R/S configuration.
    """

    # --- Problem Definition ---
    # The question asks for the EXO product of the Diels-Alder reaction between
    # 2,5-dimethylthiophene and Furan-2,5-dione (maleic anhydride).
    
    # The options as provided in the question text.
    options = {
        'A': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        'B': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        'C': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        'D': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione"
    }

    # The final answer provided by the LLM to be checked.
    llm_answer = 'A'

    # --- Verification Constraints ---

    # Constraint 1: Structural Correctness
    # The diene is 2,5-dimethylthiophene, which contains sulfur.
    # Therefore, the resulting bicyclic adduct must have a sulfur bridge, named "epithio".
    # An "epoxy" bridge (oxygen) would be incorrect.
    correct_bridge_type = "epithio"
    
    # The dienophile is maleic anhydride, which forms an "isobenzofuran-1,3-dione" core.
    correct_base_name = "isobenzofuran-1,3-dione"

    # Constraint 2: Stereochemical Correctness
    # The question asks for the EXO product. Heat favors the thermodynamically stable EXO product.
    # Through Cahn-Ingold-Prelog (CIP) analysis, the EXO isomer has the stereochemistry (3aR,4S,7R,7aS).
    # The ENDO isomer would have the stereochemistry (3aR,4R,7S,7aS).
    exo_stereochemistry = "(3aR,4S,7R,7aS)"
    endo_stereochemistry = "(3aR,4R,7S,7aS)"

    # --- Checking the LLM's Answer ---
    
    # Retrieve the full text of the chosen answer.
    answer_text = options.get(llm_answer)

    if not answer_text:
        return f"Invalid Answer. The provided answer '{llm_answer}' does not correspond to any of the options A, B, C, or D."

    # Check Constraint 1: Structural Correctness
    if correct_bridge_type not in answer_text:
        return (f"Incorrect. The answer '{llm_answer}' is wrong because the diene is thiophene, which requires a sulfur bridge named '{correct_bridge_type}'. "
                f"The chosen answer incorrectly describes an 'epoxy' (oxygen) bridge.")
    
    if correct_base_name not in answer_text:
        return (f"Incorrect. The answer '{llm_answer}' is wrong because the product's base name is incorrect. "
                f"It should be '{correct_base_name}' but is something else.")

    # Check Constraint 2: Stereochemical Correctness
    if not answer_text.startswith(exo_stereochemistry):
        # Check if it matches the ENDO stereochemistry instead.
        if answer_text.startswith(endo_stereochemistry):
            return (f"Incorrect. The answer '{llm_answer}' is wrong because the question asks for the EXO product. "
                    f"The chosen answer has the stereochemistry '{endo_stereochemistry}', which corresponds to the ENDO product, not the EXO product ('{exo_stereochemistry}').")
        else:
            # This case is unlikely given the options but is good for robustness.
            return f"Incorrect. The stereochemistry of answer '{llm_answer}' does not match the expected EXO or ENDO products."

    # If all constraints are satisfied, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_answer_correctness()
print(result)