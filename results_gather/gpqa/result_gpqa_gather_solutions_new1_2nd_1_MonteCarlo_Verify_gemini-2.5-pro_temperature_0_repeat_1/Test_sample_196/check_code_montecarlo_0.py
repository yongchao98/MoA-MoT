def check_answer():
    """
    Checks the correctness of the LLM's answer by:
    1. Identifying the starting material from spectral data.
    2. Predicting the reaction product.
    3. Comparing the predicted product with the LLM's chosen option.
    """

    # --- Data from the Question ---
    # IR: 3400–2500 (broad), 1720 -> Aromatic COOH
    # NMR: 10.5 (COOH), 8.0/7.2 (para-sub), alkyl signals
    # Alkyl NMR: 0.9(t,3H), 1.4(d,3H), 1.7(m,2H), 2.9(m,1H) -> sec-butyl
    
    # The options as presented in the final question block
    options = {
        'A': "1-isobutyl-4-methylbenzene",
        'B': "2-(4-ethylphenyl)propanoic acid",
        'C': "4-(sec-butyl)benzoic acid",
        'D': "1-(sec-butyl)-4-methylbenzene"
    }
    
    llm_final_answer_key = 'D'

    # --- Step 1: Identify the starting material (Compound X) ---
    # Based on the analysis, the starting material must be an aromatic carboxylic acid
    # with a para-substituted sec-butyl group.
    expected_starting_material = "4-(sec-butyl)benzoic acid"
    
    # Check if this starting material is among the options
    found_starting_material = False
    for key, name in options.items():
        if name == expected_starting_material:
            found_starting_material = True
            break
            
    if not found_starting_material:
        return f"Constraint Failure: The correct starting material, '{expected_starting_material}', identified from the spectral data, is not present in the options list."

    # --- Step 2: Analyze the reaction ---
    # Reagents: Red Phosphorus (P) and Hydroiodic Acid (HI)
    # Transformation: This is a strong reduction that converts a carboxylic acid (-COOH)
    # to a methyl group (-CH3).
    
    # --- Step 3: Determine the final product ---
    # Applying the transformation to the starting material:
    # 4-(sec-butyl)benzoic acid -> 1-(sec-butyl)-4-methylbenzene
    correct_product_name = "1-(sec-butyl)-4-methylbenzene"

    # --- Step 4: Check the LLM's answer ---
    llm_chosen_product_name = options.get(llm_final_answer_key)

    if not llm_chosen_product_name:
        return f"Invalid Answer: The key '{llm_final_answer_key}' provided by the LLM does not correspond to any of the options A, B, C, or D."

    if llm_chosen_product_name == correct_product_name:
        return "Correct"
    else:
        reason = (
            f"The provided answer '{llm_final_answer_key}' is incorrect.\n"
            f"1. **Starting Material Identification**: The IR and NMR data correctly identify the starting material as 4-(sec-butyl)benzoic acid (Option C).\n"
            f"2. **Reaction Prediction**: The reaction with red P and HI reduces the carboxylic acid group (-COOH) to a methyl group (-CH3).\n"
            f"3. **Correct Final Product**: The final product is therefore 1-(sec-butyl)-4-methylbenzene (Option D).\n"
            f"4. **Error**: The LLM chose option {llm_final_answer_key} ('{llm_chosen_product_name}'), which is not the correct product."
        )
        return reason

# Execute the check
result = check_answer()
print(result)