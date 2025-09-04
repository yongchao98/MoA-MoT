def check_answer():
    """
    Checks the correctness of the LLM's answer by deducing the chemical formula from the given spectroscopic data.
    """
    # --- Step 1: Define constraints from the question ---
    is_di_substituted_aromatic = True
    has_ester_group = True
    aromatic_h_signals = 2
    vinyl_h_signals = 2
    ch3_signals = 2
    ch2_signals = 0
    
    # The options provided in the question
    options = {
        "A": "C11H12O2",
        "B": "C11H14O2",
        "C": "C12H12O2",
        "D": "C12H14O2"
    }
    
    # The answer provided by the LLM
    llm_answer_option = "A"

    # --- Step 2: Deduce the structure and formula from constraints ---
    
    # Start with the base aromatic ring. Di-substituted implies C6H4.
    # 2 aromatic H signals is consistent with para-substitution.
    c_atoms = 6
    h_atoms = 4
    o_atoms = 0
    
    # Add the ester group (-COO-). This adds C and O2.
    if has_ester_group:
        c_atoms += 1
        o_atoms += 2
    else:
        return "Incorrect. The question states an ester group is present, which was not considered."

    # Identify substituents from NMR data.
    # 2 vinyl H signals + 1 CH3 signal strongly suggests a propenyl group (-CH=CH-CH3).
    # This group has the formula C3H5.
    # The constraint of no -CH2 groups is met.
    has_propenyl_group = (vinyl_h_signals == 2 and ch3_signals >= 1 and ch2_signals == 0)
    if has_propenyl_group:
        c_atoms += 3  # C3 for the propenyl group
        h_atoms += 5  # H5 for the propenyl group
        ch3_signals -= 1 # Account for the CH3 used in the propenyl group
    else:
        return "Incorrect. The logic fails to identify the propenyl group from the vinyl and CH3 signals."

    # The remaining signal is one CH3 group.
    # This must be the alkyl part of the ester, forming a methyl ester (-COOCH3).
    # The carbonyl carbon is already counted. We just add the methyl part.
    if ch3_signals == 1:
        c_atoms += 1  # C for the methyl group of the ester
        h_atoms += 3  # H3 for the methyl group of the ester
    else:
        return f"Incorrect. After accounting for the propenyl group, there should be 1 CH3 signal left, but the calculation resulted in {ch3_signals}."

    calculated_formula = f"C{c_atoms}H{h_atoms}O{o_atoms}"

    # --- Step 3: Compare the deduced formula with the LLM's answer ---
    
    # Check if the calculated formula matches any of the options.
    if calculated_formula not in options.values():
        return f"Logic error: The derived formula {calculated_formula} does not match any of the given options."

    # Check if the LLM's chosen option corresponds to the correct formula.
    correct_option_formula = options.get(llm_answer_option)
    
    if correct_option_formula == calculated_formula:
        return "Correct"
    else:
        return (f"Incorrect. The correct formula based on the spectral data is {calculated_formula}. "
                f"The LLM chose option {llm_answer_option}, which corresponds to {correct_option_formula}.")

# Run the check
result = check_answer()
print(result)