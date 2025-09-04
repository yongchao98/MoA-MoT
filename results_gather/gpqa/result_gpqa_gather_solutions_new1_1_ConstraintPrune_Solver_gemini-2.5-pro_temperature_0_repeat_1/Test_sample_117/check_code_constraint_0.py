def check_chemistry_answer():
    """
    This function checks the correctness of the answer to a chemistry question
    by encoding the rules of the reaction.
    """
    # 1. Define the problem statement and the provided answer.
    question = {
        "reactant": "4,4-dimethylcyclopent-1-enol",
        "reagent": "bromine",
        "options": {
            "A": "4-bromo-4,4-dimethylcyclopentanone",
            "B": "(1R,2R)-1,2-dibromo-4,4-dimethylcyclopentanol",
            "C": "2-bromo-4,4-dimethylcyclopentanone",
            "D": "(1R,2S)-1,2-dibromo-4,4-dimethylcyclopentanol"
        }
    }
    llm_answer_key = "C"

    # 2. Apply chemical principles to determine the correct product.

    # Principle 1: Identify the functional group and reaction type.
    # The reactant "4,4-dimethylcyclopent-1-enol" is an enol.
    # The reaction of an enol with a halogen (bromine) is a classic alpha-halogenation.
    # The major pathway is not simple addition, but substitution at the alpha-carbon
    # with rearrangement to the more stable keto form.
    major_pathway = "alpha-halogenation"
    expected_product_class = "alpha-halo-ketone"

    # Principle 2: Determine the specific structure of the product.
    # Reactant: 4,4-dimethylcyclopent-1-enol
    # - The double bond is between C1 and C2.
    # - The hydroxyl group (-OH) is on C1.
    # - The alpha-carbon (adjacent to the C-OH) is C2.
    # Product formation:
    # - The -OH at C1 becomes a ketone (C=O).
    # - The bromine atom adds to the alpha-carbon, C2.
    # - The 4,4-dimethyl groups are unchanged.
    # Correct product name: "2-bromo-4,4-dimethylcyclopentanone"
    correct_product_name = "2-bromo-4,4-dimethylcyclopentanone"

    # 3. Map the correct product to the given options.
    correct_option_key = None
    for key, value in question["options"].items():
        if value == correct_product_name:
            correct_option_key = key
            break
    
    if not correct_option_key:
        return "Error in checker: The chemically correct product was not found in the options."

    # 4. Validate the provided answer.
    if llm_answer_key == correct_option_key:
        # The answer is correct. Let's double-check why the others are wrong.
        errors = []
        if "4-bromo" in question["options"]["A"]:
            pass # Correctly identified as wrong position
        else: errors.append("Logic for A is flawed.")
        
        if "dibromo" in question["options"]["B"] and "dibromo" in question["options"]["D"]:
            pass # Correctly identified as wrong product class (addition vs substitution)
        else: errors.append("Logic for B/D is flawed.")

        if not errors:
            return "Correct"
        else:
            return f"Answer is correct, but checker logic has issues: {errors}"
    else:
        return (f"Incorrect. The provided answer is '{llm_answer_key}', but the correct answer is '{correct_option_key}'.\n"
                f"Reason: The reaction of an enol with bromine is an alpha-halogenation, which yields an alpha-halo-ketone. "
                f"The correct product is '{correct_product_name}', which corresponds to option {correct_option_key}.")

# Execute the check and print the result.
result = check_chemistry_answer()
print(result)