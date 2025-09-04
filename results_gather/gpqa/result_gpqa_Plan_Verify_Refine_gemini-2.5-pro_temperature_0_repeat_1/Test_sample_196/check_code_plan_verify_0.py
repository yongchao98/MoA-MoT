import re

def check_answer():
    """
    This function checks the correctness of the LLM's answer by verifying its reasoning process.
    1. It checks if the spectral data correctly identifies the starting material (Compound X).
    2. It verifies the chemical transformation caused by the specified reagents.
    3. It confirms that the final product matches the LLM's selected option.
    """
    
    # --- Data from the question ---
    ir_data = {
        "O-H_stretch_acid": (2500, 3400), # Broad
        "C=O_stretch_acid": 1720,
        "Aromatic_stretch": [1610, 1450]
    }
    
    nmr_data = {
        "COOH": {"shift": 10.5, "integration": 1, "multiplicity": "bs"},
        "Aromatic_ortho_COOH": {"shift": 8.0, "integration": 2, "multiplicity": "d"},
        "Aromatic_meta_COOH": {"shift": 7.2, "integration": 2, "multiplicity": "d"},
        "Benzylic_CH": {"shift": 2.9, "integration": 1, "multiplicity": "m"},
        "Alkyl_CH2": {"shift": 1.7, "integration": 2, "multiplicity": "m"},
        "Alkyl_CH3_d": {"shift": 1.4, "integration": 3, "multiplicity": "d"},
        "Alkyl_CH3_t": {"shift": 0.9, "integration": 3, "multiplicity": "t"}
    }
    
    reagent = "red phosphorus and HI"
    llm_answer_choice = "D"
    options = {
        "A": "4-(sec-butyl)benzoic acid",
        "B": "2-(4-ethylphenyl)propanoic acid",
        "C": "1-isobutyl-4-methylbenzene",
        "D": "1-(sec-butyl)-4-methylbenzene"
    }

    # --- Step 1: Verify the structure of Compound X from spectral data ---
    
    # Check IR for carboxylic acid and aromatic ring
    if not (ir_data["O-H_stretch_acid"] and ir_data["C=O_stretch_acid"]):
        return "Incorrect analysis: The IR data strongly suggests a carboxylic acid, which the reasoning should identify."
    if not ir_data["Aromatic_stretch"]:
        return "Incorrect analysis: The IR data suggests an aromatic ring, which the reasoning should identify."

    # Check NMR for key functional groups
    if not (nmr_data["COOH"]["shift"] > 10 and nmr_data["COOH"]["integration"] == 1):
        return "Incorrect analysis: The NMR signal at 10.5 ppm (1H) is characteristic of a carboxylic acid proton, a key feature of Compound X."
        
    # Check for para-disubstituted benzene ring
    if not (nmr_data["Aromatic_ortho_COOH"]["multiplicity"] == 'd' and nmr_data["Aromatic_meta_COOH"]["multiplicity"] == 'd' and \
            nmr_data["Aromatic_ortho_COOH"]["integration"] == 2 and nmr_data["Aromatic_meta_COOH"]["integration"] == 2):
        return "Incorrect analysis: The two doublets in the aromatic region (8.0 ppm and 7.2 ppm, each 2H) are a classic pattern for a 1,4- (para) disubstituted benzene ring."

    # Check for the sec-butyl group structure from the remaining NMR signals
    # -CH(a)-CH3(b)
    #   |
    #   CH2(c)-CH3(d)
    # (a) Benzylic CH: multiplet, coupled to CH3(b) (3H) and CH2(c) (2H). Total neighbors = 5. Correct.
    # (b) CH3: doublet, coupled to CH(a) (1H). Correct. (Signal at 1.4 ppm)
    # (c) CH2: multiplet, coupled to CH(a) (1H) and CH3(d) (3H). Total neighbors = 4. Correct.
    # (d) CH3: triplet, coupled to CH2(c) (2H). Correct. (Signal at 0.9 ppm)
    
    # Verify the assignment of signals to the sec-butyl group
    if not (nmr_data["Benzylic_CH"]["integration"] == 1 and nmr_data["Alkyl_CH2"]["integration"] == 2 and \
            nmr_data["Alkyl_CH3_d"]["integration"] == 3 and nmr_data["Alkyl_CH3_t"]["integration"] == 3):
        return "Incorrect analysis: The integration of the alkyl signals (1H, 2H, 3H, 3H) does not match a sec-butyl group."
    
    # If all checks pass, Compound X is correctly identified.
    compound_x_name = "4-(sec-butyl)benzoic acid"
    if options["A"] != compound_x_name:
        return f"Constraint check failed: Option A should be the starting material, {compound_x_name}, but it is listed as {options['A']}."

    # --- Step 2: Verify the reaction ---
    # Red P / HI is a strong reducing agent that reduces carboxylic acids to alkanes.
    # -COOH -> -CH3
    
    reactant_functional_group = "-COOH"
    product_functional_group = "-CH3"
    
    if "benzoic acid" in compound_x_name and reagent == "red phosphorus and HI":
        # The "benzoic acid" part becomes "methylbenzene" (or toluene)
        expected_product_name = compound_x_name.replace("benzoic acid", "methylbenzene")
        # Standardize the name: 4-(sec-butyl)methylbenzene is 1-(sec-butyl)-4-methylbenzene
        expected_product_name = "1-(sec-butyl)-4-methylbenzene"
    else:
        return "Incorrect reaction prediction: The reaction of a carboxylic acid with red P / HI should result in complete reduction to an alkane (-CH3)."

    # --- Step 3: Verify the final answer ---
    llm_answer_name = options.get(llm_answer_choice)
    
    if llm_answer_name is None:
        return f"Invalid option: The LLM chose '{llm_answer_choice}', which is not a valid option."

    if llm_answer_name == expected_product_name:
        # Final check on other options
        if options["A"] == compound_x_name: # A is starting material
            if "isobutyl" in options["C"]: # C has wrong isomer
                return "Correct"
            else:
                return "Logic check failed: The reasoning for why option C is incorrect is flawed."
        else:
            return "Logic check failed: The reasoning for why option A is incorrect is flawed."
    else:
        return f"Incorrect final product: The predicted product is '{expected_product_name}', but the LLM's answer is '{llm_answer_name}' (Option {llm_answer_choice})."

# Run the check
result = check_answer()
if result == "Correct":
    print("Correct")
else:
    print(f"The answer is incorrect.\nReason: {result}")
