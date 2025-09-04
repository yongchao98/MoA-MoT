import re

def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer by verifying each step of its reasoning.
    1. Verifies that the spectral data matches the identified starting material (Compound X).
    2. Verifies the chemical reaction logic (reduction of a carboxylic acid).
    3. Verifies that the final product matches the reaction outcome and the selected option.
    """
    # --- Define the problem and the LLM's answer ---
    # The question's data
    ir_data = ["3400–2500", "1720", "1610", "1450"]
    nmr_data = [
        {"ppm": 10.5, "multiplicity": "bs", "integration": 1}, # bs = broad singlet
        {"ppm": 8.0, "multiplicity": "d", "integration": 2},
        {"ppm": 7.2, "multiplicity": "d", "integration": 2},
        {"ppm": 2.9, "multiplicity": "m", "integration": 1},
        {"ppm": 1.7, "multiplicity": "m", "integration": 2},
        {"ppm": 1.4, "multiplicity": "d", "integration": 3},
        {"ppm": 0.9, "multiplicity": "t", "integration": 3},
    ]
    options = {
        "A": "1-(sec-butyl)-4-methylbenzene",
        "B": "4-(sec-butyl)benzoic acid",
        "C": "1-isobutyl-4-methylbenzene",
        "D": "2-(4-ethylphenyl)propanoic acid"
    }
    llm_final_choice = "A"

    # The LLM's reasoning steps
    llm_identified_compound_x = "4-(sec-butyl)benzoic acid"
    llm_predicted_reaction = "Red P / HI reduces a carboxylic acid (-COOH) to a methyl group (-CH3)"
    llm_derived_product = "1-(sec-butyl)-4-methylbenzene"

    # --- Step 1: Check the identification of the starting material (Compound X) ---

    # 1a. Check IR data consistency
    # The broad 3400-2500 and 1720 peaks indicate a carboxylic acid.
    # The 1610 and 1450 peaks indicate an aromatic ring.
    # The identified compound X, 4-(sec-butyl)benzoic acid, is an aromatic carboxylic acid, which is consistent.
    if "benzoic acid" not in llm_identified_compound_x:
        return "Incorrect: The identified starting material, '{}', is not an aromatic carboxylic acid, which contradicts the IR data (3400-2500 cm-1 and 1720 cm-1 peaks).".format(llm_identified_compound_x)

    # 1b. Check NMR data consistency
    # We will check if the NMR signals (integration and multiplicity) match the proposed structure.
    
    # Expected NMR pattern for 4-(sec-butyl)benzoic acid:
    # - COOH: 1H, singlet
    # - Aromatic: 4H total, as two 2H doublets (para-substitution)
    # - sec-butyl group [-CH(CH3)(CH2CH3)]:
    #   - Benzylic CH: 1H, multiplet (coupled to 3H+2H)
    #   - CH3 on CH: 3H, doublet (coupled to 1H)
    #   - CH2 of ethyl: 2H, multiplet (coupled to 1H+3H)
    #   - CH3 of ethyl: 3H, triplet (coupled to 2H)
    expected_nmr_counts = {
        "s_1H": 1, "d_2H": 2, "m_1H": 1, "m_2H": 1, "d_3H": 1, "t_3H": 1
    }

    # Parse the given NMR data into the same format
    given_nmr_counts = {}
    for signal in nmr_data:
        # Normalize multiplicity: 'bs' (broad singlet) is treated as 's'
        mult = signal["multiplicity"].replace('b', '')
        integ = signal["integration"]
        key = f"{mult}_{integ}H"
        given_nmr_counts[key] = given_nmr_counts.get(key, 0) + 1

    if given_nmr_counts != expected_nmr_counts:
        return "Incorrect: The NMR data does not match the proposed starting material, '{}'. Expected signal counts were {}, but the data provides {}.".format(llm_identified_compound_x, expected_nmr_counts, given_nmr_counts)

    # --- Step 2: Check the reaction prediction ---
    
    # The reaction is with Red Phosphorus and HI.
    # The LLM correctly states this is a strong reducing agent that converts -COOH to -CH3.
    # Let's apply this transformation to the identified starting material.
    reactant_name = llm_identified_compound_x # "4-(sec-butyl)benzoic acid"
    # The "benzoic acid" part implies a benzene ring with a COOH group. Reducing it makes it "toluene" (methylbenzene).
    # So, "4-(sec-butyl)benzoic acid" becomes "4-(sec-butyl)toluene".
    # The IUPAC name for 4-(sec-butyl)toluene is 1-(sec-butyl)-4-methylbenzene.
    correct_product_name = "1-(sec-butyl)-4-methylbenzene"

    if llm_derived_product != correct_product_name:
        return "Incorrect: The reaction logic is flawed. Reducing '{}' should yield '{}', but the answer derived '{}'.".format(reactant_name, correct_product_name, llm_derived_product)

    # --- Step 3: Check the final product identification ---

    # The LLM correctly derived the product name. Now check if it matches the chosen option.
    if options.get(llm_final_choice) != llm_derived_product:
        return "Incorrect: The final derived product is '{}', but the selected option {} corresponds to '{}'.".format(llm_derived_product, llm_final_choice, options.get(llm_final_choice))

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer_correctness()
print(result)