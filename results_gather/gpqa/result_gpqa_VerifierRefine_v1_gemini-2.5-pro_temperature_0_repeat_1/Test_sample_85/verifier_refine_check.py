def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer by codifying the chemical principles.
    """
    errors = []

    # --- Principle 1: Reagent Selectivity ---
    # LiBH4 selectively reduces esters over carboxylic acids.
    # BH3 selectively reduces carboxylic acids over esters.
    # The LLM's explanation correctly states this. We'll codify this for our logic.
    reagent_selectivity = {
        'LiBH4': {'reduces': 'ester', 'spares': 'carboxylic_acid'},
        'BH3': {'reduces': 'carboxylic_acid', 'spares': 'ester'}
    }

    # --- Principle 2: Stereochemical Outcome ---
    # The starting material is 3-ethyl-5-isobutoxy-5-oxopentanoic acid.
    # The chiral center is at C3 (the carbon with the ethyl group).
    # The reactions (reduction and lactonization) happen at C1 (acid) and C5 (ester).
    # Since the reactions do not involve the chiral center, its configuration is retained.
    stereochemistry_rule = "retained"

    # --- Problem Definition ---
    # Reaction A: A + LiBH4 -> (R)-product
    # Reaction B: B + BH3 -> (S)-product
    problem_definition = {
        'A': {'reagent': 'LiBH4', 'product_config': 'R'},
        'B': {'reagent': 'BH3', 'product_config': 'S'}
    }

    # --- LLM's Answer to be Checked ---
    # The LLM concluded:
    # A = (R)-3-ethyl-5-isobutoxy-5-oxopentanoic acid
    # B = (S)-3-ethyl-5-isobutoxy-5-oxopentanoic acid
    # This corresponds to option A.
    llm_answer = {
        'A_config': 'R',
        'B_config': 'S'
    }
    llm_selected_option = 'A'

    # --- Verification Logic ---

    # 1. Deduce the required configuration for starting material A
    # Reaction A produces (R)-product.
    # The stereochemistry is retained.
    # Therefore, starting material A must be (R).
    required_A_config = problem_definition['A']['product_config']
    if llm_answer['A_config'] != required_A_config:
        errors.append(
            f"Error in logic for starting material A. "
            f"The reaction produces an (R)-product and retains stereochemistry, so the starting material A must be (R). "
            f"The answer incorrectly identifies it as ({llm_answer['A_config']})."
        )

    # 2. Deduce the required configuration for starting material B
    # Reaction B produces (S)-product.
    # The stereochemistry is retained.
    # Therefore, starting material B must be (S).
    required_B_config = problem_definition['B']['product_config']
    if llm_answer['B_config'] != required_B_config:
        errors.append(
            f"Error in logic for starting material B. "
            f"The reaction produces an (S)-product and retains stereochemistry, so the starting material B must be (S). "
            f"The answer incorrectly identifies it as ({llm_answer['B_config']})."
        )
        
    # 3. Check if the final selected option matches the derived logic
    options = {
        "A": {'A_config': 'R', 'B_config': 'S'},
        "B": {'A_config': 'S', 'B_config': 'S'},
        "C": {'A_config': 'R', 'B_config': 'R'},
        "D": {'A_config': 'S', 'B_config': 'R'}
    }
    
    correctly_derived_answer = {'A_config': required_A_config, 'B_config': required_B_config}
    
    if options.get(llm_selected_option) != correctly_derived_answer:
        errors.append(
            f"The final selected option '{llm_selected_option}' is inconsistent with the correct chemical reasoning. "
            f"The correct reasoning implies A should be ({required_A_config}) and B should be ({required_B_config}), "
            f"which corresponds to option A, not {llm_selected_option}."
        )


    # --- Final Verdict ---
    if not errors:
        return "Correct"
    else:
        return "Incorrect. Reason(s):\n" + "\n".join(errors)

# Execute the check and print the result
result = check_chemistry_answer()
print(result)