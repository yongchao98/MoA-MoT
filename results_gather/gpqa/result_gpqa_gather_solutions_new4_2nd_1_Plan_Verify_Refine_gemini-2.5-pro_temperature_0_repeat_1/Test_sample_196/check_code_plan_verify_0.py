import re

def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer by simulating the chemical problem.
    1. Identifies the starting material from spectral data.
    2. Determines the reaction product.
    3. Compares the derived product with the given answer.
    """
    # --- Problem Definition ---
    # Spectral data clues for Compound X
    ir_features = {
        "carboxylic_acid_oh": "3400-2500 cm-1 (broad)",
        "carbonyl_cooh": "1720 cm-1"
    }
    nmr_features = {
        "cooh_proton": {"ppm": 10.5, "integrates": 1},
        "para_doublets": [{"ppm": 8.0, "integrates": 2}, {"ppm": 7.2, "integrates": 2}],
        "alkyl_protons": {"total_H": 1+2+3+3} # Sum of integrations for alkyl part
    }
    # Reaction details
    reagents = "red phosphorus and HI"
    
    # Provided options from the question
    options = {
        'A': "2-(4-ethylphenyl)propanoic acid",
        'B': "4-(sec-butyl)benzoic acid",
        'C': "1-(sec-butyl)-4-methylbenzene",
        'D': "1-isobutyl-4-methylbenzene"
    }
    
    # The final answer provided by the LLM to be checked
    llm_final_answer = "C"

    # --- Step 1: Identify the Starting Material (Compound X) ---
    analysis_log = []
    
    # Check for carboxylic acid
    has_carboxylic_acid = "carbonyl_cooh" in ir_features and "cooh_proton" in nmr_features
    if not has_carboxylic_acid:
        return "Constraint Failure: The spectral data clearly indicates a carboxylic acid, but this was not identified."
    analysis_log.append("Identified functional group: Carboxylic Acid (-COOH).")

    # Check for para-disubstituted benzene ring
    is_para_substituted = len(nmr_features["para_doublets"]) == 2
    if not is_para_substituted:
        return "Constraint Failure: The two doublets in the aromatic region (8.0, 7.2 ppm) indicate a 1,4-disubstituted (para) ring."
    analysis_log.append("Identified core structure: 1,4-disubstituted benzene ring.")

    # Check for sec-butyl group based on proton count and structure clues
    # A butyl group has formula C4H9, so 9 alkyl protons.
    is_butyl_group = nmr_features["alkyl_protons"]["total_H"] == 9
    # The specific splitting pattern (doublet, triplet, multiplets) points to sec-butyl over other isomers like isobutyl.
    alkyl_group = "sec-butyl" if is_butyl_group else "unknown"
    if alkyl_group != "sec-butyl":
        return "Constraint Failure: The alkyl signals in the NMR correspond to a sec-butyl group, not any other isomer."
    analysis_log.append(f"Identified alkyl substituent: {alkyl_group}.")

    # Assemble the starting material name
    starting_material_name = f"4-({alkyl_group})benzoic acid"
    analysis_log.append(f"Deduced starting material (Compound X): {starting_material_name}.")

    # Verify that the deduced starting material is one of the options (a common distractor)
    if starting_material_name != options['B']:
        return f"Logic Error: The deduced starting material '{starting_material_name}' does not match option B, which it should."

    # --- Step 2: Determine the Reaction Product ---
    product_name = ""
    if reagents == "red phosphorus and HI":
        # This is a known, powerful reduction of a carboxylic acid to a methyl group.
        # The transformation is from "benzoic acid" to "methylbenzene" (toluene).
        if "benzoic acid" in starting_material_name:
            # Replace the "benzoic acid" part with "methylbenzene" and adjust numbering
            base_name = starting_material_name.replace("benzoic acid", "methylbenzene")
            # The name changes from 4-(...) to 1-(...)-4-...
            product_name = base_name.replace("4-(", "1-(").replace(")", ")-4")
            analysis_log.append(f"Reaction: Red P/HI reduces -COOH to -CH3.")
            analysis_log.append(f"Deduced final product: {product_name}.")
        else:
            return "Reaction Logic Error: Expected a benzoic acid derivative as starting material."
    else:
        return f"Unknown Reagents: The logic for reagents '{reagents}' is not defined."

    # --- Step 3: Match Product to Options and Verify Answer ---
    correct_option_letter = None
    for letter, name in options.items():
        if name == product_name:
            correct_option_letter = letter
            break
    
    if correct_option_letter is None:
        return f"Final Check Failure: The derived product '{product_name}' was not found in the provided options."

    if llm_final_answer == correct_option_letter:
        return "Correct"
    else:
        return (f"Incorrect. The final answer is wrong.\n"
                f"Reason: The analysis correctly identifies the starting material as '{starting_material_name}' (Option B) "
                f"and the final product as '{product_name}'. This product corresponds to Option {correct_option_letter}, "
                f"but the provided answer was Option {llm_final_answer}.")

# Run the check
result = check_chemistry_answer()
print(result)