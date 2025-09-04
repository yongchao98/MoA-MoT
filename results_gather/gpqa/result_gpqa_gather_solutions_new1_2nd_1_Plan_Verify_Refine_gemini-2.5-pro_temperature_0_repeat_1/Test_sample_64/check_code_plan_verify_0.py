def check_chemistry_answer():
    """
    Checks the correctness of the proposed answer based on chemical principles.
    """
    # Define the multiple-choice options from the question
    options = {
        "A": {
            "A": "decahydro-7H-benzo[7]annulen-7-one",
            "B": "3-ethylpent-4-enoic acid"
        },
        "B": {
            "A": "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one",
            "B": "lithium 3-ethylpent-4-enoate"
        },
        "C": {
            "A": "decahydro-7H-benzo[7]annulen-7-one",
            "B": "lithium 3-ethylpent-4-enoate"
        },
        "D": {
            "A": "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one",
            "B": "3-ethylpent-4-enoic acid"
        }
    }

    # The final answer to be checked
    proposed_answer_key = "B"

    # Known correct information based on chemical principles
    correct_info = {
        "start_material_A_formula": "C11H16O",
        "product_A_formulas": {
            "decahydro-7H-benzo[7]annulen-7-one": "C11H18O",
            "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one": "C11H16O"
        },
        "correct_product_B_form": "lithium 3-ethylpent-4-enoate" # The salt form
    }

    # Retrieve the specific answer choice to check
    answer_to_check = options.get(proposed_answer_key)

    if not answer_to_check:
        return f"Invalid answer key '{proposed_answer_key}'. Key must be one of {list(options.keys())}."

    # Constraint 1: Check molecular formula for product A
    # The reaction is an isomerization, so the product formula must match the starting material.
    product_A_name = answer_to_check["A"]
    product_A_formula = correct_info["product_A_formulas"].get(product_A_name)
    
    if product_A_formula != correct_info["start_material_A_formula"]:
        return (f"Incorrect. Constraint on Product A is not satisfied. "
                f"The reaction for A is an isomerization, so the product's molecular formula must be "
                f"{correct_info['start_material_A_formula']}. The proposed product A, '{product_A_name}', "
                f"has a formula of {product_A_formula}, which is incorrect.")

    # Constraint 2: Check the form of product B
    # The reaction uses a strong lithium base (LDA) and has no specified acidic workup,
    # so the product must be the lithium salt, not the free acid.
    product_B_name = answer_to_check["B"]
    if product_B_name != correct_info["correct_product_B_form"]:
        return (f"Incorrect. Constraint on Product B is not satisfied. "
                f"The reaction for B has no specified acidic workup. Therefore, the product should be the "
                f"lithium salt ('{correct_info['correct_product_B_form']}'), not the free acid ('{product_B_name}').")

    # If both constraints are satisfied, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_chemistry_answer()
print(result)