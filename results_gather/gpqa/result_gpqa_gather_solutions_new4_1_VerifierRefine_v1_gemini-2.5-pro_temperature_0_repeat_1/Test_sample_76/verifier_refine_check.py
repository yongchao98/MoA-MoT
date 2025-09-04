def check_chemistry_answer():
    """
    Checks the correctness of the selected answer based on chemical principles.
    """
    
    # The final answer provided by the LLM to be checked.
    llm_answer = 'D'

    # Define the products for each option as given in the question.
    options = {
        'A': {
            'A': "(Z)-2-methyl-5-phenylpent-2-en-1-ol",
            'B': "2,3,4,6,7,8-hexamethyl-5,9,10,11-tetrahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorine"
        },
        'B': {
            'A': "4-methyl-1-phenylpent-3-en-1-ol",
            'B': "2,3,4,6,7,8-hexamethyl-5,9,10,11-tetrahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorine"
        },
        'C': {
            'A': "(Z)-2-methyl-5-phenylpent-2-en-1-ol",
            'B': "2,3,4,6,7,8-hexamethyl-5,9,9a,10,11,11a-hexahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorene"
        },
        'D': {
            'A': "4-methyl-1-phenylpent-3-en-1-ol",
            'B': "2,3,4,6,7,8-hexamethyl-5,9,9a,10,11,11a-hexahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorene"
        }
    }

    # --- Step 1: Retrieve the products for the selected answer ---
    selected_products = options.get(llm_answer)
    if not selected_products:
        return f"Invalid answer choice '{llm_answer}'. The choice must be one of {list(options.keys())}."

    product_A_name = selected_products['A']
    product_B_name = selected_products['B']

    # --- Step 2: Check the product of Reaction A ---
    # The [1,2]-Wittig rearrangement is the most plausible pathway, leading to '4-methyl-1-phenylpent-3-en-1-ol'.
    # The other option is structurally inconsistent with this rearrangement.
    correct_product_A_name = "4-methyl-1-phenylpent-3-en-1-ol"
    if product_A_name != correct_product_A_name:
        return (f"Incorrect: Product A is wrong. The Wittig rearrangement should yield "
                f"'{correct_product_A_name}', but the answer proposes '{product_A_name}'.")

    # --- Step 3: Check the product of Reaction B ---
    # The Cope rearrangement is an isomerization, so the degree of saturation must be conserved.
    # The starting material is a 'hexahydro' derivative, so the product must also be 'hexahydro'.
    if 'hexahydro' not in product_B_name:
        # Determine if the proposed product is 'tetrahydro' or something else for a more specific error message.
        saturation_level = 'tetrahydro' if 'tetrahydro' in product_B_name else 'different'
        return (f"Incorrect: Product B is wrong. The Cope rearrangement is an isomerization, so the 'hexahydro' "
                f"starting material must yield a 'hexahydro' product. The answer proposes a '{saturation_level}' product.")

    # --- Step 4: Final Conclusion ---
    # If both checks pass, the answer is correct.
    return "Correct"

# Run the checker and print the result
result = check_chemistry_answer()
print(result)