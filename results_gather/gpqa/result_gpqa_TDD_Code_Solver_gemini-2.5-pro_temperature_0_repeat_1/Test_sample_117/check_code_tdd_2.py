def check_organic_reaction_answer():
    """
    This function checks the correctness of the answer to a specific organic chemistry question.
    It uses a knowledge-based approach to determine the expected product.
    """
    question = "What is the major outcome of the reaction between 4,4-dimethylcyclopent-1-enol and bromine?"
    llm_answer = "<<<B>>>"

    # --- Knowledge Base and Analysis ---

    # 1. Define reactants and reaction type
    reactant_name = "4,4-dimethylcyclopent-1-enol"
    reagent_name = "bromine"

    # 2. Analyze the reactant structure and reaction type
    # The name "enol" indicates the presence of a C=C double bond with an -OH group on one of the carbons.
    # The reaction of an enol with a halogen (Br2) is a classic alpha-halogenation.
    # The enol tautomerizes to its keto form, and the halogen adds to the alpha-carbon.
    
    # Structure breakdown of 4,4-dimethylcyclopent-1-enol:
    # - cyclopent: 5-membered ring
    # - -1-en-1-ol: The double bond is between C1 and C2, and the hydroxyl is on C1.
    # - 4,4-dimethyl: Two methyl groups on C4.

    # Prediction of the reaction outcome:
    # - The enol group (C1=C2-OH) will become a ketone at C1 (C1=O).
    # - The bromine will add to the alpha-carbon, which was C2 of the enol double bond.
    # - The rest of the structure (4,4-dimethyl groups) remains unchanged.
    
    # 3. Construct the expected product name
    # Base name: 4,4-dimethylcyclopentanone
    # Substituent: bromo at position 2
    expected_product_name = "2-bromo-4,4-dimethylcyclopentanone"

    # 4. Define the given options
    options = {
        "A": "(1R,2R)-1,2-dibromo-4,4-dimethylcyclopentanol",
        "B": "2-bromo-4,4-dimethylcyclopentanone",
        "C": "4-bromo-4,4-dimethylcyclopentanone",
        "D": "(1R,2S)-1,2-dibromo-4,4-dimethylcyclopentanol"
    }

    # 5. Determine the correct option key
    correct_option_key = None
    for key, value in options.items():
        if value == expected_product_name:
            correct_option_key = key
            break
    
    if correct_option_key is None:
        # This case should not be reached if the options are well-defined.
        return "Error in checker: The predicted correct product is not among the options."

    # 6. Validate the provided LLM answer
    # Extract the letter from the answer format, e.g., "<<<B>>>" -> "B"
    try:
        provided_option_key = llm_answer.strip().replace("<", "").replace(">", "")
    except:
        return f"Invalid answer format: {llm_answer}"

    if provided_option_key == correct_option_key:
        return "Correct"
    else:
        # Generate a reason for the incorrectness
        chosen_product = options.get(provided_option_key, "an invalid option")
        correct_product = options[correct_option_key]

        if provided_option_key in ["A", "D"]:
            reason = (f"The answer {provided_option_key} is incorrect. The product '{chosen_product}' would result from "
                      f"a standard electrophilic addition of Br2 to an alkene, forming a dihaloalcohol. "
                      f"However, the starting material is an enol. The reaction of an enol with bromine is an "
                      f"alpha-halogenation, which yields an alpha-haloketone. The correct product is '{correct_product}' (Option {correct_option_key}).")
        elif provided_option_key == "C":
            reason = (f"The answer {provided_option_key} is incorrect. The product '{chosen_product}' suggests bromination at the C4 position. "
                      f"Alpha-halogenation of an enol specifically occurs at the alpha-carbon (C2 in this case), not at other saturated positions on the ring. "
                      f"The correct product is '{correct_product}' (Option {correct_option_key}).")
        else:
            reason = (f"The provided answer {provided_option_key} is incorrect. The correct answer is {correct_option_key}. "
                      f"The reaction is an alpha-halogenation of an enol, leading to the product '{correct_product}'.")
        
        return reason

# Execute the check and print the result
result = check_organic_reaction_answer()
print(result)