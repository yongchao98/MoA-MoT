import re

def check_reaction_outcome(selected_answer_key: str):
    """
    Checks the correctness of the proposed product for the reaction between
    4,4-dimethylcyclopent-1-enol and bromine.

    Args:
        selected_answer_key: The key ('A', 'B', 'C', or 'D') of the proposed answer.

    Returns:
        A string "Correct" if the answer is correct, or a string explaining the error.
    """
    options = {
        "A": "2-bromo-4,4-dimethylcyclopentanone",
        "B": "4-bromo-4,4-dimethylcyclopentanone",
        "C": "(1R,2R)-1,2-dibromo-4,4-dimethylcyclopentanol",
        "D": "(1R,2S)-1,2-dibromo-4,4-dimethylcyclopentanol"
    }

    if selected_answer_key not in options:
        return f"Error: The provided answer key '{selected_answer_key}' is not a valid option."

    product_name = options[selected_answer_key]

    # Constraint 1: Keto-Enol Tautomerism.
    # The starting material, an enol, exists in equilibrium with its more stable keto tautomer.
    # The major reaction pathway is the alpha-halogenation of this ketone.
    # Therefore, the major product must be a ketone (name ends in "-one").
    if not product_name.endswith("one"):
        return (f"Incorrect. The product '{product_name}' is an alcohol (functional group '-ol'). "
                "The starting enol tautomerizes to the more stable ketone, 4,4-dimethylcyclopentanone. "
                "The major reaction is the alpha-halogenation of this ketone, which must result in a ketone product.")

    # Constraint 2: Regioselectivity of Alpha-Halogenation.
    # The reaction is an alpha-halogenation, meaning the bromine atom replaces a hydrogen
    # on a carbon adjacent to the carbonyl group (C=O).
    # In 4,4-dimethylcyclopentanone, the carbonyl is at C1.
    # The alpha-carbons (adjacent carbons) are C2 and C5.
    # The product must have bromine at one of these positions.
    
    # Use regex to find the position number of the bromo group.
    match = re.search(r'(\d+)-bromo', product_name)
    
    if not match:
        return f"Error: Could not parse the bromine position from the name '{product_name}'."
        
    bromo_position = int(match.group(1))
    
    alpha_positions = [2, 5]
    if bromo_position not in alpha_positions:
        return (f"Incorrect. The product '{product_name}' has bromine at position {bromo_position}. "
                "This is an alpha-halogenation reaction, so the bromine must substitute at an alpha-position. "
                f"For 4,4-dimethylcyclopentanone, the alpha-positions relative to the C1 carbonyl are {alpha_positions}.")

    # If both constraints are satisfied, the answer is correct.
    return "Correct"

# The LLM's answer is "A". Let's run the check.
llm_answer = "A"
result = check_reaction_outcome(llm_answer)
print(result)

# You can also test why other options are wrong:
# print("\nChecking other options for verification:")
# print(f"B: {check_reaction_outcome('B')}")
# print(f"C: {check_reaction_outcome('C')}")
# print(f"D: {check_reaction_outcome('D')}")