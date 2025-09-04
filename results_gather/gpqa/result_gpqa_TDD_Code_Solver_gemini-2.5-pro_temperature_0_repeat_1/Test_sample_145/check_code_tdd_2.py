def check_diels_alder_stereochemistry():
    """
    Checks the correctness of the selected answer for the Diels-Alder reaction
    between 5-fluorocyclopenta-1,3-diene and maleic anhydride.

    The code verifies the answer against two fundamental stereochemical constraints:
    1.  Endo Selectivity: The major product must be the 'endo' adduct. For this
        specific molecular structure, 'endo' corresponds to the '(4S,7R)'
        stereochemical descriptor.
    2.  Syn Selectivity: Due to the electronic 'syn-directing' effect of the
        electronegative fluorine atom, the major product must be the 'syn'
        adduct. This corresponds to the '8s' stereochemical descriptor.

    The correct answer must satisfy both constraints.
    """
    # The answer provided by the LLM
    llm_answer = "A"

    # The options from the question
    options = {
        "A": "(3aR,4S,7R,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "B": "(3aR,4R,7S,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "C": "(3aR,4R,7S,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "D": "(3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
    }

    # Constraint 1: The product must be 'endo'.
    # This corresponds to the '(...,4S,7R,...)' descriptor.
    endo_descriptor = "(4S,7R,"

    # Constraint 2: The product must be 'syn'.
    # This corresponds to the '(...,8s)' descriptor.
    syn_descriptor = ",8s)"

    # Find the option that satisfies both constraints
    major_product_key = None
    for key, name in options.items():
        is_endo = endo_descriptor in name
        is_syn = syn_descriptor in name
        if is_endo and is_syn:
            major_product_key = key
            break

    # Check if the LLM's answer matches the chemically correct product
    if llm_answer == major_product_key:
        return "Correct"
    elif major_product_key is None:
        return "Error in checking logic: No option satisfies both the endo and syn selectivity rules."
    else:
        # Analyze the incorrectness of the LLM's choice
        llm_choice_name = options.get(llm_answer)
        if not llm_choice_name:
            return f"Incorrect. The provided answer '{llm_answer}' is not a valid option."

        is_llm_choice_endo = endo_descriptor in llm_choice_name
        is_llm_choice_syn = syn_descriptor in llm_choice_name

        reasons = []
        if not is_llm_choice_endo:
            reasons.append("it is an 'exo' adduct, but the kinetically favored product is 'endo'")
        if not is_llm_choice_syn:
            reasons.append("it is an 'anti' adduct, but the electronically favored product is 'syn'")

        return (f"Incorrect. The answer '{llm_answer}' is wrong because {', and '.join(reasons)}. "
                f"The major product must be the endo-syn adduct, which is option '{major_product_key}'.")

# Execute the check
result = check_diels_alder_stereochemistry()
print(result)