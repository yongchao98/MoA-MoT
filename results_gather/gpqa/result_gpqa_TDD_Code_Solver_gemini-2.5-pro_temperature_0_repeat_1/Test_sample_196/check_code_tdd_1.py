def check_chemistry_answer():
    """
    This function checks the correctness of the given LLM's answer to the chemistry problem.
    It follows the logical steps of spectral analysis and reaction prediction.
    """
    # --- Problem Data ---
    ir_data = ["3400–2500 cm-1", "1720 cm-1", "1610 cm-1", "1450 cm-1"]
    nmr_data = ["10.5 ppm (bs, 1H)", "8.0 ppm (d, 2H)", "7.2 ppm (d, 2H)", 
                "2.9 ppm (m, 1H)", "1.7 ppm (m, 2H)", "1.4 ppm (d, 3H)", "0.9 ppm (t, 3H)"]
    reagents = "Red P / HI"
    options = {
        "A": "4-(sec-butyl)benzoic acid",
        "B": "2-(4-ethylphenyl)propanoic acid",
        "C": "1-isobutyl-4-methylbenzene",
        "D": "1-(sec-butyl)-4-methylbenzene"
    }
    llm_answer_key = "D"

    # --- Step 1: Verify the identity of the starting material (Compound X) ---
    
    # Check IR for carboxylic acid features
    is_carboxylic_acid_ir = "3400–2500 cm-1" in ir_data and "1720 cm-1" in ir_data
    if not is_carboxylic_acid_ir:
        return "Incorrect: The IR data strongly suggests a carboxylic acid, but this feature was not correctly identified."

    # Check NMR for key features of 4-(sec-butyl)benzoic acid
    nmr_signals = set(nmr_data)
    expected_nmr_features = {
        "cooh_proton": "10.5 ppm (bs, 1H)",
        "para_aromatic_1": "8.0 ppm (d, 2H)",
        "para_aromatic_2": "7.2 ppm (d, 2H)",
        "sec_butyl_ch": "2.9 ppm (m, 1H)",
        "sec_butyl_ch2": "1.7 ppm (m, 2H)",
        "sec_butyl_ch3_d": "1.4 ppm (d, 3H)",
        "sec_butyl_ch3_t": "0.9 ppm (t, 3H)"
    }

    for feature, signal in expected_nmr_features.items():
        if signal not in nmr_signals:
            return f"Incorrect: The starting material identification is flawed. The NMR signal for '{feature}' ({signal}) is missing from the provided data."

    identified_starting_material = "4-(sec-butyl)benzoic acid"
    if identified_starting_material != options["A"]:
        return f"Incorrect: The identified starting material '{identified_starting_material}' should correspond to option A, but it does not."

    # --- Step 2: Predict the product of the reaction ---
    
    # Red P / HI reduces a carboxylic acid (-COOH) to a methyl group (-CH3)
    if reagents == "Red P / HI" and "benzoic acid" in identified_starting_material:
        # The name "4-(sec-butyl)benzoic acid" transforms to "1-(sec-butyl)-4-methylbenzene"
        predicted_product = "1-(sec-butyl)-4-methylbenzene"
    else:
        return "Incorrect: The reaction prediction logic is flawed. Red P/HI should reduce the carboxylic acid."

    # --- Step 3: Compare with the LLM's answer ---
    
    # Find the option key for the predicted product
    correct_key = None
    for key, value in options.items():
        if value == predicted_product:
            correct_key = key
            break
    
    if correct_key is None:
        return f"Incorrect: The correctly predicted product '{predicted_product}' is not found in the options."

    if llm_answer_key == correct_key:
        return "Correct"
    else:
        return (f"Incorrect: The LLM's answer is {llm_answer_key}, but the correct answer is {correct_key}. "
                f"The starting material is {identified_starting_material}, which upon reduction with {reagents} "
                f"yields {predicted_product}.")

result = check_chemistry_answer()
print(result)