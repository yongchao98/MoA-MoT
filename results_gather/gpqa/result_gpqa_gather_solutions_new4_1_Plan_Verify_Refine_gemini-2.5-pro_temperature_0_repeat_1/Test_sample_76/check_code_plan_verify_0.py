def check_chemistry_answer():
    """
    Checks the correctness of the selected option based on chemical principles.
    """
    # The final answer provided by the LLM to be checked.
    llm_answer = 'D'

    # Define the products for each option from the question.
    options = {
        'A': {
            'A_product_name': "(Z)-2-methyl-5-phenylpent-2-en-1-ol",
            'B_product_name': "2,3,4,6,7,8-hexamethyl-5,9,10,11-tetrahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorine"
        },
        'B': {
            'A_product_name': "4-methyl-1-phenylpent-3-en-1-ol",
            'B_product_name': "2,3,4,6,7,8-hexamethyl-5,9,10,11-tetrahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorine"
        },
        'C': {
            'A_product_name': "(Z)-2-methyl-5-phenylpent-2-en-1-ol",
            'B_product_name': "2,3,4,6,7,8-hexamethyl-5,9,9a,10,11,11a-hexahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorene"
        },
        'D': {
            'A_product_name': "4-methyl-1-phenylpent-3-en-1-ol",
            'B_product_name': "2,3,4,6,7,8-hexamethyl-5,9,9a,10,11,11a-hexahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorene"
        }
    }

    # --- Define Correctness Criteria ---

    # Criterion 1: Product A must be the result of the most plausible Wittig rearrangement pathway.
    # This is the [1,2]-rearrangement product.
    correct_A_name = "4-methyl-1-phenylpent-3-en-1-ol"

    # Criterion 2: Product B must be an isomer of the starting material.
    # The Cope rearrangement conserves the molecular formula and degree of saturation.
    # The starting material is 'hexahydro', so the product must also be 'hexahydro'.
    required_saturation_in_B = "hexahydro"

    # --- Check the LLM's Answer ---
    selected_products = options.get(llm_answer)

    if not selected_products:
        return f"Invalid option '{llm_answer}' was provided. Please choose from A, B, C, or D."

    # Check if the selected option satisfies Criterion 1 for Reaction A
    if selected_products['A_product_name'] != correct_A_name:
        return (f"Incorrect: The product for Reaction A is wrong. "
                f"The most plausible pathway for the Wittig rearrangement leads to '{correct_A_name}', "
                f"but the selected answer chose '{selected_products['A_product_name']}'.")

    # Check if the selected option satisfies Criterion 2 for Reaction B
    if required_saturation_in_B not in selected_products['B_product_name']:
        return (f"Incorrect: The product for Reaction B is wrong. "
                f"The Cope rearrangement is an isomerization, so the 'hexahydro' starting material must yield a 'hexahydro' product. "
                f"The selected answer's product for B is not described as 'hexahydro'.")

    # If all criteria are met, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_chemistry_answer()
print(result)