import re

def check_correctness():
    """
    This function checks the correctness of the provided answer to a chemistry question.
    It codifies the rules of the reaction between an enol and a halogen.
    """
    # --- Problem Definition ---
    question = "What is the major outcome of the reaction between 4,4-dimethylcyclopent-1-enol and bromine?"
    options = {
        "A": "(1R,2S)-1,2-dibromo-4,4-dimethylcyclopentanol",
        "B": "2-bromo-4,4-dimethylcyclopentanone",
        "C": "4-bromo-4,4-dimethylcyclopentanone",
        "D": "(1R,2R)-1,2-dibromo-4,4-dimethylcyclopentanol"
    }
    # The final answer provided by the LLM analysis.
    llm_final_answer = "B"

    # --- Chemical Rule Engine ---

    # Rule 1: Determine the major reaction pathway.
    # The reaction of an enol with a halogen (like Br2) is a classic alpha-substitution.
    # The driving force is the formation of the thermodynamically stable carbonyl (C=O) group.
    # This means the major product must be a ketone, not an alcohol (which would be an addition product).
    
    product_name = options.get(llm_final_answer)
    if not product_name:
        return f"Invalid answer option '{llm_final_answer}'. Please choose from A, B, C, D."

    if "anol" in product_name:
        return (f"Incorrect. The answer '{product_name}' is an alcohol. "
                "This would be the result of an electrophilic addition reaction. "
                "However, the major pathway for the reaction of an enol with bromine is alpha-substitution, "
                "which leads to the formation of a more stable ketone.")

    if "ketone" not in product_name:
        return f"Incorrect. The major product should be a ketone, but the answer '{product_name}' does not appear to be one."

    # Rule 2: Determine the position of halogenation.
    # The reaction is an *alpha*-halogenation. The halogen adds to the alpha-carbon.
    # Reactant: 4,4-dimethylcyclopent-1-enol
    # The enol is at C1(-OH)=C2. The alpha-carbon is C2.
    # The corresponding ketone is 4,4-dimethylcyclopentanone (carbonyl at C1).
    # The alpha-carbons are C2 and C5, which are equivalent.
    # Therefore, bromination must occur at position 2 (or 5).
    
    expected_product = "2-bromo-4,4-dimethylcyclopentanone"

    if product_name == expected_product:
        return "Correct"
    
    # Check for common incorrect substitution patterns
    if "4-bromo" in product_name:
        return (f"Incorrect. The answer '{product_name}' suggests bromination at the C4 position. "
                "The reaction is an alpha-halogenation, which occurs at the carbon adjacent to the carbonyl group "
                "(the alpha-position, C2), not at the saturated C4 position.")

    # Generic failure message if it's a ketone but not the right one
    return f"Incorrect. While the answer '{product_name}' is a ketone, it does not match the expected product '{expected_product}'."

# Execute the check and print the result
result = check_correctness()
print(result)