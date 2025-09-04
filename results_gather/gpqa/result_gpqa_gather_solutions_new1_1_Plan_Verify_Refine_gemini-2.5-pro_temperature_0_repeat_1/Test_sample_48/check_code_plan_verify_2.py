def get_molecular_formula(compound_name):
    """
    Returns the molecular formula for a given compound name based on its structure.
    This function acts as a knowledge base for the checker.
    """
    formulas = {
        # --- Reactants ---
        "1,1-dimethoxyethan-1-amine": "C4H11NO2",
        "but-3-en-2-ol": "C4H8O",
        "(3R,4S)-3,4-dimethylhexa-1,5-diyne": "C8H10",
        "2-((vinyloxy)methyl)but-1-ene": "C7H12O",
        
        # --- Products for A ---
        "(Z)-1-(but-2-en-2-yloxy)ethen-1-amine": "C6H11NO", # CH2=C(NH2)-O-C(CH3)=CH-CH3
        "6-methyl-3,4-dihydro-2H-pyran-2-amine": "C6H11NO",
        
        # --- Products for B ---
        "(3Z,4E)-3,4-diethylidenecyclobut-1-ene": "C8H10", # C4H2(ring) + 2*C2H4(ethylidene)
        "(1Z,2E)-1,2-diethylidenecyclobutane": "C8H12", # C4H6(ring) + 2*C2H4(ethylidene)
        
        # --- Products for C ---
        "4-methylenehexanal": "C7H12O", # O=CH-CH2-CH2-C(=CH2)-CH2-CH3
        "4-methylenehexan-1-ol": "C7H14O", # HO-CH2-CH2-CH2-C(=CH2)-CH2-CH3
    }
    return formulas.get(compound_name, "Unknown")

def check_answer():
    """
    Checks the correctness of the final answer by verifying chemical constraints.
    """
    llm_answer = 'A' # The answer provided by the last LLM to be checked.

    options = {
        'A': {'A': "(Z)-1-(but-2-en-2-yloxy)ethen-1-amine", 'B': "(3Z,4E)-3,4-diethylidenecyclobut-1-ene", 'C': "4-methylenehexanal"},
        'B': {'A': "6-methyl-3,4-dihydro-2H-pyran-2-amine", 'B': "(1Z,2E)-1,2-diethylidenecyclobutane", 'C': "4-methylenehexan-1-ol"},
        'C': {'A': "6-methyl-3,4-dihydro-2H-pyran-2-amine", 'B': "(1Z,2E)-1,2-diethylidenecyclobutane", 'C': "4-methylenehexanal"},
        'D': {'A': "(Z)-1-(but-2-en-2-yloxy)ethen-1-amine", 'B': "(3Z,4E)-3,4-diethylidenecyclobut-1-ene", 'C': "4-methylenehexan-1-ol"}
    }

    # --- Constraint 1: Reaction B is an isomerization ---
    # The product must have the same formula as the reactant.
    reactant_b_formula = get_molecular_formula("(3R,4S)-3,4-dimethylhexa-1,5-diyne")
    product_b_formula_A = get_molecular_formula(options['A']['B'])
    product_b_formula_B = get_molecular_formula(options['B']['B'])
    
    if product_b_formula_A != reactant_b_formula:
        return f"Incorrect. The provided answer is '{llm_answer}', but option A is invalid. Reason: Product B in option A, '{options['A']['B']}', has formula {product_b_formula_A}, which does not match the reactant formula {reactant_b_formula}."
    if product_b_formula_B == reactant_b_formula:
        return f"Incorrect. The provided answer is '{llm_answer}', but this is wrong. Reason: The product for reaction B must be an isomer of the reactant (C8H10). The product in options B and C, '{options['B']['B']}', has formula {product_b_formula_B} (C8H12), which is incorrect. The product in options A and D, '{options['A']['B']}', has formula {product_b_formula_A} (C8H10), which is correct. This eliminates options B and C."

    # --- Constraint 2: Reaction C is an isomerization ---
    # The product must have the same formula as the reactant.
    reactant_c_formula = get_molecular_formula("2-((vinyloxy)methyl)but-1-ene")
    product_c_formula_A = get_molecular_formula(options['A']['C'])
    product_c_formula_D = get_molecular_formula(options['D']['C'])

    if product_c_formula_A != reactant_c_formula:
        return f"Incorrect. The provided answer is '{llm_answer}', but option A is invalid. Reason: Product C in option A, '{options['A']['C']}', has formula {product_c_formula_A}, which does not match the reactant formula {reactant_c_formula}."
    if product_c_formula_D == reactant_c_formula:
        return f"Incorrect. The provided answer is '{llm_answer}', but this is wrong. Reason: The product for reaction C must be an isomer of the reactant (C7H12O). The product in options A and C, '{options['A']['C']}', has formula {product_c_formula_A} (C7H12O), which is correct. The product in options B and D, '{options['D']['C']}', has formula {product_c_formula_D} (C7H14O), which is incorrect. This eliminates options B and D."

    # --- Constraint 3: Reaction C product type ---
    # A Claisen rearrangement produces a carbonyl, not an alcohol.
    if "ol" in options['D']['C']:
        pass # This is already covered by the formula check, but it's a good chemical check.
    else:
        return f"Incorrect. The provided answer is '{llm_answer}', but this is wrong. Reason: The product for reaction C, a Claisen rearrangement, must be a carbonyl compound (aldehyde/ketone). The product in options B and D is an alcohol, which is incorrect."

    # --- Final Verification ---
    # Based on the checks above:
    # - Reaction B eliminates options B and C.
    # - Reaction C eliminates options B and D.
    # The only option that satisfies all constraints is A.
    
    if llm_answer == 'A':
        return "Correct"
    else:
        return f"Incorrect. The provided answer is '{llm_answer}'. Based on atom conservation constraints, the only valid option is A. Reaction B eliminates options B and C (product formula C8H12 vs reactant C8H10). Reaction C eliminates options B and D (product formula C7H14O vs reactant C7H12O)."

# Run the check
result = check_answer()
print(result)