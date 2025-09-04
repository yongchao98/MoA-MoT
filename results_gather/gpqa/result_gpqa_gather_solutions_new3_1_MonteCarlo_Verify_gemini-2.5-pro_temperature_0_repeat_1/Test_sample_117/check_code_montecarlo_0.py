import re

def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer to a chemistry question.

    The function encodes the chemical principles for the reaction of an enol with a halogen
    to determine the correct product and compares it with the provided answer.
    """
    # --- Problem Definition ---
    question = "What is the major outcome of the reaction between 4,4-dimethylcyclopent-1-enol and bromine?"
    options = {
        "A": "(1R,2R)-1,2-dibromo-4,4-dimethylcyclopentanol",
        "B": "(1R,2S)-1,2-dibromo-4,4-dimethylcyclopentanol",
        "C": "4-bromo-4,4-dimethylcyclopentanone",
        "D": "2-bromo-4,4-dimethylcyclopentanone"
    }
    llm_final_answer = "<<<D>>>"

    # --- Chemical Analysis ---
    # 1. Identify reactants and reaction type.
    # Reactant: 4,4-dimethylcyclopent-1-enol. The key functional group is an enol (C=C-OH).
    # Reagent: Bromine (Br2). A halogen that acts as an electrophile.
    # Reaction: Reaction of an enol with a halogen.

    # 2. Determine the major reaction pathway.
    # For an enol reacting with a halogen, the major pathway is alpha-halogenation.
    # This is because the reaction is driven by the formation of the very stable
    # carbonyl (C=O) bond in the final product. Electrophilic addition to form a
    # dihaloalcohol is a minor pathway.

    # 3. Predict the product of alpha-halogenation.
    # - The starting enol is 4,4-dimethylcyclopent-1-enol. The double bond is between C1 and C2,
    #   and the -OH group is on C1.
    # - The "alpha-carbon" is the carbon of the double bond that does NOT bear the -OH group. This is C2.
    # - The bromine atom adds to the alpha-carbon (C2).
    # - The enol group at C1 tautomerizes to its more stable keto form (a C=O group).
    # - The 4,4-dimethyl groups remain unchanged.
    # - The resulting structure is a cyclopentanone ring, with a carbonyl at C1, a bromine at C2,
    #   and two methyl groups at C4.
    # - The name of this product is 2-bromo-4,4-dimethylcyclopentanone.

    correct_product_name = "2-bromo-4,4-dimethylcyclopentanone"

    # 4. Match the predicted product to the given options.
    correct_option_key = None
    for key, value in options.items():
        # We can ignore stereochemical descriptors like (1R,2R) for this comparison
        # as the core structure is what matters.
        normalized_option_name = re.sub(r'^\(.*\)-', '', value).strip()
        if normalized_option_name == correct_product_name:
            correct_option_key = key
            break

    if correct_option_key is None:
        return f"Checker Error: The predicted correct product '{correct_product_name}' was not found in the options."

    # 5. Extract the LLM's chosen option from its response.
    match = re.search(r'<<<([A-D])>>>', llm_final_answer)
    if not match:
        return f"Invalid answer format: The answer '{llm_final_answer}' does not follow the '<<<X>>>' format."
    
    llm_chosen_option = match.group(1)

    # 6. Compare the LLM's choice with the correct option and return the result.
    if llm_chosen_option == correct_option_key:
        return "Correct"
    else:
        reason = (
            f"The provided answer '{llm_chosen_option}' is incorrect.\n"
            f"The correct answer is '{correct_option_key}'.\n"
            f"Reasoning: The reaction of an enol (4,4-dimethylcyclopent-1-enol) with bromine (Br2) is an alpha-halogenation. "
            f"The major product is formed by adding a bromine atom to the alpha-carbon (C2) and converting the enol at C1 into a ketone. "
            f"This yields 2-bromo-4,4-dimethylcyclopentanone, which corresponds to option {correct_option_key}. "
            f"The answer chose option {llm_chosen_option} ('{options[llm_chosen_option]}'), which is incorrect."
        )
        return reason

# Execute the check and print the result.
result = check_answer_correctness()
print(result)