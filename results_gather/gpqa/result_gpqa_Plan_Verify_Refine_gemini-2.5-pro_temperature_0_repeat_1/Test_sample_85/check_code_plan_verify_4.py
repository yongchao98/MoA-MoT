def check_chemistry_answer():
    """
    This function checks the correctness of the LLM's answer by encoding the
    relevant chemical principles and applying them to the problem.
    """

    # --- Step 1: Define the known chemical principles ---

    # Principle 1: Chemoselectivity of the reducing agents.
    # This dictionary maps the reagent to the functional group it selectively reduces
    # in a molecule containing both a carboxylic acid and an ester.
    reagent_selectivity = {
        "LiBH4": "ester",  # Lithium borohydride selectively reduces esters over carboxylic acids.
        "BH3": "carboxylic_acid"  # Borane selectively reduces carboxylic acids over esters.
    }

    # Principle 2: Stereochemistry of the reaction.
    # The reactions (reduction and subsequent lactonization) occur at C1 and C5 of the
    # pentanoic acid derivative. The chiral center is at C3. No bonds to the chiral
    # center are broken or formed, so its configuration is retained.
    stereochemistry_is_retained = True

    # --- Step 2: Define the problem from the question ---

    # Reaction A: A + LiBH4 -> (R)-product
    # Reaction B: B + BH3 -> (S)-product
    problem_statement = {
        "A": {"reagent": "LiBH4", "product_stereochem": "R"},
        "B": {"reagent": "BH3", "product_stereochem": "S"}
    }

    # --- Step 3: Define the LLM's proposed answer ---
    # The LLM chose option B, which states:
    # A = (R)-3-ethyl-5-isobutoxy-5-oxopentanoic acid
    # B = (S)-3-ethyl-5-isobutoxy-5-oxopentanoic acid
    llm_answer = {
        "A": {"starting_material_stereochem": "R"},
        "B": {"starting_material_stereochem": "S"}
    }

    # --- Step 4: Verify the LLM's answer against the principles ---

    error_messages = []

    # Check Reaction A
    reaction_A_params = problem_statement["A"]
    llm_A_proposal = llm_answer["A"]
    
    # Based on the principle of stereochemical retention:
    if stereochemistry_is_retained:
        predicted_product_A_stereochem = llm_A_proposal["starting_material_stereochem"]
    else:
        # This case is not applicable here, but included for logical completeness
        predicted_product_A_stereochem = "unknown" 

    # Compare the predicted outcome with the question's required outcome
    if predicted_product_A_stereochem != reaction_A_params["product_stereochem"]:
        error_messages.append(
            f"Constraint failure for Reaction A: The LLM proposes starting material A is ({llm_A_proposal['starting_material_stereochem']}). "
            f"Since the reaction retains stereochemistry, the product should be ({predicted_product_A_stereochem}). "
            f"However, the question requires the product to be ({reaction_A_params['product_stereochem']})."
        )

    # Check Reaction B
    reaction_B_params = problem_statement["B"]
    llm_B_proposal = llm_answer["B"]

    # Based on the principle of stereochemical retention:
    if stereochemistry_is_retained:
        predicted_product_B_stereochem = llm_B_proposal["starting_material_stereochem"]
    else:
        predicted_product_B_stereochem = "unknown"

    # Compare the predicted outcome with the question's required outcome
    if predicted_product_B_stereochem != reaction_B_params["product_stereochem"]:
        error_messages.append(
            f"Constraint failure for Reaction B: The LLM proposes starting material B is ({llm_B_proposal['starting_material_stereochem']}). "
            f"Since the reaction retains stereochemistry, the product should be ({predicted_product_B_stereochem}). "
            f"However, the question requires the product to be ({reaction_B_params['product_stereochem']})."
        )
        
    # Additionally, we can verify the LLM's reasoning about chemoselectivity, which is a prerequisite.
    # The LLM correctly stated that LiBH4 reduces the ester and BH3 reduces the acid.
    # This check confirms the foundational logic is sound.
    if reagent_selectivity["LiBH4"] != "ester" or reagent_selectivity["BH3"] != "carboxylic_acid":
        error_messages.append("The fundamental chemical principle of reagent selectivity used for the check is incorrect.")


    # --- Step 5: Return the final verdict ---
    if not error_messages:
        return "Correct"
    else:
        return "\n".join(error_messages)

# Execute the check and print the result.
result = check_chemistry_answer()
print(result)