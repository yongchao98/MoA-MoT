def check_chemistry_answer():
    """
    Checks the correctness of the answer to the organic chemistry question.
    """
    # --- Problem Definition ---
    question = "What is the major outcome of the reaction between 4,4-dimethylcyclopent-1-enol and bromine?"
    options = {
        'A': '4-bromo-4,4-dimethylcyclopentanone',
        'B': '(1R,2S)-1,2-dibromo-4,4-dimethylcyclopentanol',
        'C': '2-bromo-4,4-dimethylcyclopentanone',
        'D': '(1R,2R)-1,2-dibromo-4,4-dimethylcyclopentanol'
    }
    llm_answer_key = 'C'

    # --- Chemical Analysis Logic ---
    
    # 1. Analyze the starting material and reagent
    reactant = "4,4-dimethylcyclopent-1-enol"
    reagent = "bromine"

    # 2. Determine the reaction type based on functional groups
    # The reactant is an enol (C=C-OH). The reagent is a halogen.
    # The reaction of an enol with a halogen is a classic alpha-halogenation.
    # The major pathway is alpha-substitution, not simple addition, due to the
    # formation of the highly stable carbonyl (C=O) group.
    
    # 3. Predict the structure of the major product
    # Reactant: 4,4-dimethylcyclopent-1-enol
    #   - The enol is at C1(-OH) and C2.
    #   - The corresponding keto tautomer is 4,4-dimethylcyclopentanone (carbonyl at C1).
    #   - The alpha-carbons (next to the carbonyl) are C2 and C5.
    # Reaction: Alpha-halogenation with bromine.
    #   - A bromine atom substitutes a hydrogen on an alpha-carbon.
    #   - The positions C2 and C5 are equivalent due to symmetry.
    # Product: Bromination at C2 gives 2-bromo-4,4-dimethylcyclopentanone.
    expected_product_name = "2-bromo-4,4-dimethylcyclopentanone"

    # --- Verification ---
    
    # Check if the key for the expected product exists in the options
    correct_key = None
    for key, value in options.items():
        if value == expected_product_name:
            correct_key = key
            break
    
    if correct_key is None:
        return f"Error in question setup: The expected correct product '{expected_product_name}' is not among the options."

    # Check if the LLM's chosen answer is correct
    if llm_answer_key == correct_key:
        # Further check if the other options are indeed incorrect for the right reasons.
        # Check A: Bromination at C4 is incorrect.
        if "4-bromo" in options['A']:
            pass # Correctly identified as a distractor.
        else:
            return "Logic Error: Option A is not the expected distractor."

        # Check B and D: Dibromo-alcohols are products of addition, not the major substitution pathway.
        if "dibromo" in options['B'] and "dibromo" in options['D']:
            pass # Correctly identified as distractors.
        else:
            return "Logic Error: Options B/D are not the expected distractors."
            
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is '{options[llm_answer_key]}' (Option {llm_answer_key}), "
                f"but the major product of the reaction is '{expected_product_name}' (Option {correct_key}). "
                f"The reaction is an alpha-halogenation of an enol, which results in an alpha-bromo ketone, not an addition product or substitution at a non-alpha position.")

# Run the check
result = check_chemistry_answer()
print(result)