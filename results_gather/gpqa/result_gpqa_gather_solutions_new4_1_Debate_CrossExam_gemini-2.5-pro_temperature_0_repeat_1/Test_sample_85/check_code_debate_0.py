def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer by applying chemical principles.
    1.  Chemoselectivity of LiBH4 and BH3.
    2.  Stereochemical outcome based on changes in Cahn-Ingold-Prelog (CIP) priorities.
    """
    
    # The final answer provided by the LLM to be checked.
    llm_answer_option = "A"

    # --- Define Chemical Principles ---

    # 1. Chemoselectivity of reducing agents
    # This dictionary defines which functional group is preferentially reduced.
    reagent_selectivity = {
        "LiBH4": "reduces_ester",
        "BH3": "reduces_acid"
    }

    # 2. Cahn-Ingold-Prelog (CIP) priorities of relevant groups.
    # A lower number means a higher priority. This is a simplified model
    # sufficient to determine the relative order.
    cip_priority = {
        "ester_side": 1,   # -CH2COOiBu
        "acid_side": 2,    # -CH2COOH
        "alcohol_side": 3, # -CH2CH2OH
        "ethyl_group": 4   # -CH2CH3
    }

    # --- Analysis of Reaction A ---
    # A + LiBH4 + H+ ---> (R)-product
    
    # Pathway: Ester is reduced to alcohol.
    # Groups on chiral center before: ester_side, acid_side, ethyl_group
    # Groups on chiral center after: alcohol_side, acid_side, ethyl_group
    
    # Check if the priority order of the top two groups is swapped.
    # Before: ester_side (1) > acid_side (2)
    # After: acid_side (2) > alcohol_side (3)
    # The group that was #1 is now #2, and the group that was #2 is now #1. This is an inversion.
    stereochemical_change_A = "inversion"
    
    target_product_A = "R"
    if stereochemical_change_A == "inversion":
        required_start_A = "S"
    else: # retention
        required_start_A = "R"

    # --- Analysis of Reaction B ---
    # B + BH3 + H+ ---> (S)-product

    # Pathway: Acid is reduced to alcohol.
    # Groups on chiral center before: ester_side, acid_side, ethyl_group
    # Groups on chiral center after: ester_side, alcohol_side, ethyl_group

    # Check if the priority order of the top two groups is maintained.
    # Before: ester_side (1) > acid_side (2)
    # After: ester_side (1) > alcohol_side (3)
    # The highest priority group remains the highest priority group. This is retention.
    stereochemical_change_B = "retention"

    target_product_B = "S"
    if stereochemical_change_B == "retention":
        required_start_B = "S"
    else: # inversion
        required_start_B = "R"

    # --- Final Verification ---
    
    # Map the LLM's answer option to the configurations
    option_map = {
        "A": {"A": "S", "B": "S"},
        "B": {"A": "R", "B": "S"},
        "C": {"A": "S", "B": "R"},
        "D": {"A": "R", "B": "R"}
    }

    if llm_answer_option not in option_map:
        return f"Error: The provided answer '{llm_answer_option}' is not a valid option (A, B, C, or D)."

    llm_config = option_map[llm_answer_option]
    errors = []

    # Compare derived requirements with the LLM's answer
    if llm_config["A"] != required_start_A:
        error_msg = (
            f"Constraint Violated for Starting Material A: "
            f"The answer states A is ({llm_config['A']}), but analysis requires A to be ({required_start_A}).\n"
            f"Reason: The reaction with LiBH4 reduces the ester. This causes the CIP priorities of the groups on the chiral center to swap, "
            f"leading to an *inversion* of the R/S descriptor. To get an (R)-product, an (S)-starting material is needed."
        )
        errors.append(error_msg)
    
    if llm_config["B"] != required_start_B:
        error_msg = (
            f"Constraint Violated for Starting Material B: "
            f"The answer states B is ({llm_config['B']}), but analysis requires B to be ({required_start_B}).\n"
            f"Reason: The reaction with BH3 reduces the carboxylic acid. This does not change the CIP priority order, "
            f"leading to a *retention* of the R/S descriptor. To get an (S)-product, an (S)-starting material is needed."
        )
        errors.append(error_msg)

    if errors:
        return "\n".join(errors)
    else:
        return "Correct"

# Execute the check
result = check_chemistry_answer()
print(result)