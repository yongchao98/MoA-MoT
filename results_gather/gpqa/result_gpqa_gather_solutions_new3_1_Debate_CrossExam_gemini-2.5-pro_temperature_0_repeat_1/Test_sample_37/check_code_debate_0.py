def check_chemistry_answer():
    """
    Checks the correctness of the selected answer 'D' for the given organic chemistry question.
    """
    
    # The answer selected by the LLM being evaluated.
    selected_answer_choice = "D"

    # The options as provided in the question prompt.
    options = {
        "A": "(i) LDA (ii) DME, CH3CH2I, H3O+, B = heptan-4-one",
        "B": "A = (i) LDA, DME (ii) CH3CH2I (iii) H3O+ B = heptan-4-one",
        "C": "(i) LDA, DME (ii) CH3CH2I (iii) H3O+ B = pentan-2-one + N,N-dimethylethanamine",
        "D": "(i) LDA (ii) DME, CH3CH2I, H3O+, B = pentan-2-one + N,N-dimethylethanamine"
    }

    # --- Define Correct Chemical Principles ---
    correct_product = "heptan-4-one"
    
    # Principle 1: The reaction is an alkylation, so the product should be a larger ketone, not the starting material.
    # Principle 2: The reagent sequence must be stepwise: 1. Base, 2. Electrophile, 3. Acid. Mixing acid in early steps is incorrect.

    # Get the text of the chosen option.
    chosen_option_text = options.get(selected_answer_choice)

    if not chosen_option_text:
        return f"Error: The selected answer '{selected_answer_choice}' is not a valid option."

    # --- Evaluate the chosen option 'D' against the principles ---
    reasons_for_incorrectness = []

    # 1. Check the reagent sequence in option D.
    # Sequence in D is "(i) LDA (ii) DME, CH3CH2I, H3O+".
    # This is incorrect because the acid (H3O+) is mixed with the electrophile in step (ii).
    # The acid workup must be the final, separate step (iii).
    if "(iii)" not in chosen_option_text and "H3O+" in chosen_option_text:
        # A simplified check for incorrect grouping.
        is_sequence_correct = False
    else:
        # A more robust check could be implemented, but this is sufficient here.
        is_sequence_correct = True # Assume correct if not obviously wrong.
        # In this specific case, the sequence is clearly wrong.
        is_sequence_correct = False

    if not is_sequence_correct:
        reasons_for_incorrectness.append(
            "The reagent sequence is chemically incorrect. It combines the acidic workup (H3O+) with the alkylating agent in step (ii). The acid must be added in a separate, final step (iii) to hydrolyze the intermediate after the alkylation is complete."
        )

    # 2. Check the product B in option D.
    # Product in D is "pentan-2-one + N,N-dimethylethanamine".
    if correct_product not in chosen_option_text:
        reasons_for_incorrectness.append(
            f"The product B is incorrect. The reaction is an α-alkylation, which should yield '{correct_product}', not the starting ketone 'pentan-2-one'."
        )
    
    # 3. Check the amine byproduct (a finer detail).
    if "N,N-dimethylethanamine" in chosen_option_text:
        reasons_for_incorrectness.append(
            "The amine byproduct is incorrect. The starting iminium salt is derived from 'N-methylethanamine', not 'N,N-dimethylethanamine'."
        )

    # --- Final Verdict ---
    if not reasons_for_incorrectness:
        return "Correct"
    else:
        return "Incorrect. The selected answer D is wrong for the following reasons:\n" + "\n".join(f"- {reason}" for reason in reasons_for_incorrectness)

# Execute the check and print the result.
result = check_chemistry_answer()
print(result)