def check_answer():
    """
    Checks the correctness of the LLM's answer based on chemical principles.
    """
    # Define the reaction products from the question
    product_A = "4-isopropylcyclohexane-1,2-dione"
    product_B = "5-methylhexane-2,3-dione"

    # --- Chemical Principle Analysis ---
    # The reaction (NaNO2, HCl, H2O) is an alpha-oxidation of a ketone.
    # It converts a -CH2- group next to a C=O group into another C=O group.
    
    # To get product A (a 1,2-dione), the starting material must be the corresponding 1-ketone.
    correct_reactant_A = "4-isopropylcyclohexan-1-one"
    
    # To get product B (a 2,3-dione), the starting material must be the corresponding 2-ketone.
    # The oxidation occurs at the alpha-methylene group (C3).
    correct_reactant_B = "5-methylhexan-2-one"

    # Define the options provided in the question
    options = {
        "A": {"A": "4-isopropyl-2-methoxycyclohexan-1-ol", "B": "5-methylhexan-2-one"},
        "B": {"A": "4-isopropylcyclohexan-1-one", "B": "5-methylhexan-2-one"},
        "C": {"A": "4-isopropyl-2-methoxycyclohexan-1-ol", "B": "5-methylhexane-2,3-diol"},
        "D": {"A": "4-isopropylcyclohexan-1-one", "B": "5-methylhexane-2,3-diol"}
    }

    # The final answer from the LLM to be checked
    llm_answer = "B"

    # --- Verification Logic ---

    # 1. Determine the truly correct option based on our chemical analysis
    actual_correct_option = None
    for option_key, reactants in options.items():
        if reactants["A"] == correct_reactant_A and reactants["B"] == correct_reactant_B:
            actual_correct_option = option_key
            break
    
    # 2. Compare the LLM's answer with the determined correct option
    if llm_answer == actual_correct_option:
        return "Correct"
    else:
        reason = f"The provided answer is '{llm_answer}', but the correct answer is '{actual_correct_option}'.\n"
        
        # Analyze why the LLM's chosen option is incorrect
        chosen_reactants = options[llm_answer]
        
        if chosen_reactants["A"] != correct_reactant_A:
            reason += f"Reason for A: The proposed starting material '{chosen_reactants['A']}' is incorrect. "
            if "ol" in chosen_reactants["A"]:
                reason += "It is an alcohol, but the reaction requires a ketone. "
            reason += f"The correct precursor for '{product_A}' is '{correct_reactant_A}'.\n"

        if chosen_reactants["B"] != correct_reactant_B:
            reason += f"Reason for B: The proposed starting material '{chosen_reactants['B']}' is incorrect. "
            if "diol" in chosen_reactants["B"]:
                reason += "It is a diol, but the reaction requires a ketone. "
            reason += f"The correct precursor for '{product_B}' is '{correct_reactant_B}'.\n"
            
        return reason.strip()

# Execute the check
result = check_answer()
print(result)