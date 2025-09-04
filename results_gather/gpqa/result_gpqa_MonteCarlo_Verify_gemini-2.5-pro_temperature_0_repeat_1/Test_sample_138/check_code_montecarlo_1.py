def check_correctness():
    """
    Checks the correctness of the given answer for the chemistry question.

    The reaction is the alpha-hydroxylation of a ketone using NaNO2/HCl/H2O
    to form an alpha-diketone. This means a -CH2- group adjacent to a
    carbonyl group is converted into another carbonyl group.
    """

    # --- Define Problem Constraints ---
    # Reaction A: ? ---> 4-isopropylcyclohexane-1,2-dione
    # Reaction B: ? ---> 5-methylhexane-2,3-dione
    llm_answer = "<<<C>>>"

    options = {
        "A": {"A": "4-isopropylcyclohexan-1-one", "B": "5-methylhexane-2,3-diol"},
        "B": {"A": "4-isopropyl-2-methoxycyclohexan-1-ol", "B": "5-methylhexane-2,3-diol"},
        "C": {"A": "4-isopropylcyclohexan-1-one", "B": "5-methylhexan-2-one"},
        "D": {"A": "4-isopropyl-2-methoxycyclohexan-1-ol", "B": "5-methylhexan-2-one"}
    }

    def check_reaction_A(start_material):
        """
        Checks if start_material yields 4-isopropylcyclohexane-1,2-dione.
        The precursor must be a ketone. The logical precursor is 4-isopropylcyclohexan-1-one,
        where the reaction occurs at the alpha-carbon C2.
        """
        # Constraint 1: Must be a ketone.
        if not start_material.endswith("-one"):
            return False, f"Starting material A ('{start_material}') is not a ketone, but the reaction requires a ketone."
        
        # Constraint 2: Must form the correct product.
        if start_material == "4-isopropylcyclohexan-1-one":
            return True, ""
        else:
            return False, f"Starting material A ('{start_material}') would not form 4-isopropylcyclohexane-1,2-dione."

    def check_reaction_B(start_material):
        """
        Checks if start_material yields 5-methylhexane-2,3-dione.
        The precursor must be a ketone. The logical precursors are 5-methylhexan-2-one
        (reaction at C3) or 5-methylhexan-3-one (reaction at C2).
        """
        # Constraint 1: Must be a ketone.
        if not start_material.endswith("-one"):
            return False, f"Starting material B ('{start_material}') is not a ketone, but the reaction requires a ketone."
        
        # Constraint 2: Must form the correct product.
        valid_precursors = ["5-methylhexan-2-one", "5-methylhexan-3-one"]
        if start_material in valid_precursors:
            return True, ""
        else:
            return False, f"Starting material B ('{start_material}') would not form 5-methylhexane-2,3-dione."

    # --- Evaluation ---
    # Extract the chosen option letter from the LLM's answer
    chosen_option_letter = llm_answer.strip('<>').upper()

    if chosen_option_letter not in options:
        return f"Invalid answer format or option '{chosen_option_letter}'. The answer must be one of {list(options.keys())}."

    chosen_materials = options[chosen_option_letter]
    
    # Check if the chosen option is correct
    is_A_correct, reason_A = check_reaction_A(chosen_materials["A"])
    if not is_A_correct:
        return f"The answer {chosen_option_letter} is incorrect. {reason_A}"

    is_B_correct, reason_B = check_reaction_B(chosen_materials["B"])
    if not is_B_correct:
        return f"The answer {chosen_option_letter} is incorrect. {reason_B}"

    # If both checks pass for the chosen option, we verify it's the *only* correct one.
    correct_options_found = []
    for option, materials in options.items():
        a_ok, _ = check_reaction_A(materials["A"])
        b_ok, _ = check_reaction_B(materials["B"])
        if a_ok and b_ok:
            correct_options_found.append(option)
    
    if len(correct_options_found) == 1 and chosen_option_letter == correct_options_found[0]:
        return "Correct"
    elif len(correct_options_found) == 0:
        return "Error in evaluation: No valid option was found based on the chemical principles."
    elif len(correct_options_found) > 1:
        return f"Error in question design: Multiple options {correct_options_found} are chemically correct."
    else: # The chosen option was incorrect, and we already have the reason.
        # This part is redundant due to earlier checks but serves as a fallback.
        return f"The answer {chosen_option_letter} is incorrect. The correct answer is {correct_options_found[0]}."

# Execute the check and print the result
result = check_correctness()
print(result)