def check_cycloaddition_answer():
    """
    Checks the correctness of the final answer for the given Diels-Alder reaction.

    The function verifies two main constraints:
    1.  Structural Skeleton: The product must have an "epithio" (sulfur) bridge,
        as the diene is thiophene.
    2.  Stereochemistry: The question asks for the EXO product, which is the
        thermodynamically favored isomer. For this reaction, the EXO product has
        the (3aR,4S,7R,7aS) stereochemical configuration.
    """
    # --- Problem Definition ---
    # The options as they are defined in the original question prompt.
    options = {
        'A': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        'B': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        'C': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        'D': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione"
    }

    # The final answer provided by the LLM to be checked.
    final_answer_choice = 'A'

    # --- Chemical Principles for Verification ---
    # 1. Structural Constraint: The diene is 2,5-dimethylthiophene (contains sulfur).
    #    Therefore, the product must have a sulfur bridge, named "epithio".
    correct_bridge_type = "epithio"
    
    # 2. Stereochemical Constraint: The question asks for the EXO product.
    #    Heat favors the thermodynamically stable EXO isomer.
    #    Based on Cahn-Ingold-Prelog rules for this structure, the EXO isomer
    #    has the (3aR,4S,7R,7aS) configuration.
    correct_exo_stereochem = "(3aR,4S,7R,7aS)"
    # The ENDO isomer would have the (3aR,4R,7S,7aS) configuration.
    correct_endo_stereochem = "(3aR,4R,7S,7aS)"

    # --- Checking the Final Answer ---
    if final_answer_choice not in options:
        return f"Invalid answer choice '{final_answer_choice}'. It must be one of {list(options.keys())}."

    selected_answer_text = options[final_answer_choice]

    # Check 1: Structural Skeleton (Bridge Type)
    if correct_bridge_type not in selected_answer_text:
        if "epoxy" in selected_answer_text:
            wrong_bridge = "epoxy"
        else:
            wrong_bridge = "an unknown type"
        return (f"Incorrect. The selected answer '{final_answer_choice}' describes a product with an '{wrong_bridge}' bridge. "
                f"Constraint not satisfied: The diene is thiophene, so the product must have a sulfur ('{correct_bridge_type}') bridge.")

    # Check 2: Stereochemistry (EXO vs. ENDO)
    if correct_exo_stereochem not in selected_answer_text:
        # Check if it's the ENDO product instead
        if correct_endo_stereochem in selected_answer_text:
            return (f"Incorrect. The selected answer '{final_answer_choice}' describes the ENDO product. "
                    f"Constraint not satisfied: The question specifically asks for the EXO product, which has the "
                    f"'{correct_exo_stereochem}' configuration.")
        else:
            # Some other incorrect stereochemistry
            return (f"Incorrect. The stereochemistry of the selected answer '{final_answer_choice}' does not match the "
                    f"required EXO product ('{correct_exo_stereochem}').")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_cycloaddition_answer()
print(result)