def check_chemistry_answer():
    """
    Checks the correctness of the selected answer for the given chemistry question
    by verifying key chemical principles.
    """
    # --- Problem Definition ---
    # The final answer from the LLM to be checked
    llm_final_answer = "B"

    # Define the options available in the question
    options = {
        "A": {
            "product_A": "decahydro-7H-benzo[7]annulen-7-one",
            "product_B": "3-ethylpent-4-enoic acid"
        },
        "B": {
            "product_A": "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one",
            "product_B": "lithium 3-ethylpent-4-enoate"
        },
        "C": {
            "product_A": "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one",
            "product_B": "3-ethylpent-4-enoic acid"
        },
        "D": {
            "product_A": "decahydro-7H-benzo[7]annulen-7-one",
            "product_B": "lithium 3-ethylpent-4-enoate"
        }
    }

    # --- Chemical Data for Verification ---
    # Molecular formulas are used to check for conservation of mass in rearrangements.
    molecular_formulas = {
        "1-vinylspiro[3.5]non-5-en-1-ol": "C11H16O",
        "decahydro-7H-benzo[7]annulen-7-one": "C11H18O",
        "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one": "C11H16O",
    }

    # Reaction conditions determine the final state of the product.
    reaction_2_conditions = "acetyl bromide (Base = LDA)"

    # --- Verification Logic ---
    selected_option = options.get(llm_final_answer)
    if not selected_option:
        return f"Error: The provided answer '{llm_final_answer}' is not a valid option (A, B, C, or D)."

    # Constraint 1: Check Reaction 1 (Anionic Oxy-Cope Rearrangement)
    # This is an isomerization, so the molecular formula must be conserved.
    reactant_A_formula = molecular_formulas["1-vinylspiro[3.5]non-5-en-1-ol"]
    product_A_name = selected_option["product_A"]
    product_A_formula = molecular_formulas.get(product_A_name)

    if reactant_A_formula != product_A_formula:
        return (f"Incorrect: Product A is wrong. The reaction is a rearrangement, so the molecular formula must be conserved. "
                f"The reactant's formula is {reactant_A_formula}, but the proposed product '{product_A_name}' has the formula {product_A_formula}.")

    # Constraint 2: Check Reaction 2 (Ireland-Claisen Rearrangement)
    # The product should be a salt because there is no acidic workup specified.
    product_B_name = selected_option["product_B"]
    
    # Check if an acidic workup (like H+ or H3O+) is mentioned for reaction 2
    acidic_workup_present = "H+" in reaction_2_conditions or "H3O+" in reaction_2_conditions

    is_acid = "oic acid" in product_B_name
    is_salt = "oate" in product_B_name

    if is_acid and not acidic_workup_present:
        return (f"Incorrect: Product B is wrong. The reaction uses a strong base (LDA) and no acidic workup is specified. "
                f"Therefore, the final product should be the deprotonated salt ('...-oate'), not the free carboxylic acid ('...-oic acid').")

    if is_salt and acidic_workup_present:
        return (f"Incorrect: Product B is wrong. An acidic workup is specified, which would protonate the salt to form the free acid.")

    # If both constraints are satisfied, the answer is correct.
    return "Correct"

# Run the check
result = check_chemistry_answer()
print(result)