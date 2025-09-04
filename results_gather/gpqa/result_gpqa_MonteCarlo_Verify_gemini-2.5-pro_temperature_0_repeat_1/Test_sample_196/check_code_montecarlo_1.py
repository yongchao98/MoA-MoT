import collections

def check_correctness_of_chemistry_answer():
    """
    This function checks the correctness of the LLM's answer to a chemistry question.
    It programmatically analyzes the spectral data, predicts the reaction product,
    and compares it to the given options and the LLM's answer.
    """
    # --- Problem Definition ---
    # This section encodes the information from the question.
    ir_data = {'broad_range': (3400, 2500), 'sharp_peak': 1720}
    nmr_data = [
        {'ppm': 10.5, 'type': 'bs', 'H': 1},  # Carboxylic acid proton
        {'ppm': 8.0, 'type': 'd', 'H': 2},   # Aromatic H ortho to COOH
        {'ppm': 7.2, 'type': 'd', 'H': 2},   # Aromatic H ortho to alkyl group
        {'ppm': 2.9, 'type': 'm', 'H': 1},   # Benzylic CH of sec-butyl
        {'ppm': 1.7, 'type': 'm', 'H': 2},   # -CH2- of sec-butyl
        {'ppm': 1.4, 'type': 'd', 'H': 3},   # -CH3 next to CH (sec-butyl)
        {'ppm': 0.9, 'type': 't', 'H': 3}    # -CH3 next to CH2 (sec-butyl)
    ]
    reagents = "red phosphorus and HI"
    llm_answer_key = 'C'
    options = {
        'A': '2-(4-ethylphenyl)propanoic acid',
        'B': '4-(sec-butyl)benzoic acid',
        'C': '1-(sec-butyl)-4-methylbenzene',
        'D': '1-isobutyl-4-methylbenzene'
    }

    # --- Step 1: Identify the Starting Material (Compound X) ---

    # 1a. Check for a carboxylic acid functional group from IR and NMR.
    has_cooh_ir = ir_data['broad_range'] == (3400, 2500) and ir_data['sharp_peak'] == 1720
    has_cooh_nmr = any(p['ppm'] > 10 and p['type'] == 'bs' for p in nmr_data)
    if not (has_cooh_ir and has_cooh_nmr):
        return "Analysis Error: The spectral data does not consistently indicate a carboxylic acid, which is the first step in identifying Compound X."

    # 1b. Check for a para-disubstituted benzene ring.
    aromatic_signals = [p for p in nmr_data if 6.5 < p['ppm'] < 8.5]
    is_para_substituted = (len(aromatic_signals) == 2 and
                           all(p['H'] == 2 for p in aromatic_signals) and
                           all(p['type'] == 'd' for p in aromatic_signals))
    if not is_para_substituted:
        return "Analysis Error: The NMR aromatic signals (two doublets, 2H each) strongly suggest a para-disubstituted ring, but the code failed to confirm this pattern."

    # 1c. Identify the alkyl group.
    # A sec-butyl group [-CH(CH3)(CH2CH3)] gives: 1H(m), 3H(d), 2H(m), 3H(t).
    alkyl_signals = [p for p in nmr_data if p['ppm'] < 4]
    signal_summary = collections.Counter([(p['type'], p['H']) for p in alkyl_signals])
    expected_sec_butyl_signals = collections.Counter([('m', 1), ('m', 2), ('d', 3), ('t', 3)])
    
    if signal_summary != expected_sec_butyl_signals:
        # An isobutyl group [-CH2CH(CH3)2] would give 2H(d), 1H(m), 6H(d).
        return "Analysis Error: The NMR signals for the alkyl group do not match a sec-butyl group. They are characteristic of a sec-butyl group, so this check points to a flaw in the logic, but the chemical reasoning is sound."

    # Conclusion for Compound X
    starting_material = "4-(sec-butyl)benzoic acid"

    # --- Step 2: Predict the Reaction Product ---
    # The reaction is with red phosphorus and HI, a strong reducing agent.
    # This reagent reduces a carboxylic acid (-COOH) to a methyl group (-CH3).
    if "benzoic acid" in starting_material and reagents == "red phosphorus and HI":
        predicted_product = "1-(sec-butyl)-4-methylbenzene"
    else:
        return "Analysis Error: The code failed to predict the outcome of the reaction."

    # --- Step 3: Find the Correct Option Key ---
    correct_key = None
    for key, name in options.items():
        if name == predicted_product:
            correct_key = key
            break
    
    if correct_key is None:
        return f"Analysis Error: The predicted product '{predicted_product}' was not found in the options."

    # --- Step 4: Final Verification ---
    if llm_answer_key == correct_key:
        return "Correct"
    else:
        # Check for common mistakes.
        if options.get(llm_answer_key) == starting_material:
            return f"Incorrect. The answer '{llm_answer_key}' ({options.get(llm_answer_key)}) is the starting material (Compound X), not the final product of the reaction with red P/HI."
        
        if "isobutyl" in options.get(llm_answer_key, ""):
             return f"Incorrect. The answer '{llm_answer_key}' suggests an isobutyl group. However, the NMR data (specifically the triplet at 0.9 ppm and doublet at 1.4 ppm) clearly indicates a sec-butyl group on the starting material, which is retained in the product."

        return (f"Incorrect. The provided answer is {llm_answer_key} ({options[llm_answer_key]}), but the correct answer is {correct_key} ({options[correct_key]}). "
                f"Reasoning: The starting material is {starting_material}. The reaction with red phosphorus and HI reduces the carboxylic acid group to a methyl group, yielding {predicted_product}.")

# Run the checker and print the result.
result = check_correctness_of_chemistry_answer()
print(result)