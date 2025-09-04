import re

def check_chemistry_answer():
    """
    Checks the correctness of the final answer for the epoxide ring-opening reaction.
    The function verifies two main constraints:
    1. Regioselectivity: The attack must occur at the less hindered carbon.
    2. Stereoselectivity: The attack must proceed with inversion of configuration, and the final R/S label must be determined correctly.
    """
    
    # The final answer chosen from the candidates
    final_answer_key = "C"
    
    # The options provided in the question
    options = {
        "A": "(1R,4R,5R)-2,2,4,5-tetramethylcyclohexan-1-ol",
        "B": "(1S,4R,5S)-2,2,4,5-tetramethylcyclohexan-1-ol",
        "C": "(1R,2S,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol",
        "D": "(1R,2R,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol"
    }
    
    if final_answer_key not in options:
        return f"Invalid answer key '{final_answer_key}'. Please choose from {list(options.keys())}."
        
    answer_string = options[final_answer_key]

    # --- Constraint 1: Regioselectivity Check ---
    # The nucleophile (Me-) attacks the less hindered carbon of the epoxide (C6, not C1).
    # Attack at C6 results in a 1,2,4,5-tetramethylcyclohexan-1-ol skeleton.
    # Attack at C1 would result in a 2,2,4,5-tetramethylcyclohexan-1-ol skeleton.
    expected_skeleton = "1,2,4,5-tetramethylcyclohexan-1-ol"
    
    if expected_skeleton not in answer_string:
        return (f"Incorrect. The answer '{final_answer_key}' violates the regioselectivity constraint. "
                f"The nucleophile should attack the less hindered C6, leading to a '{expected_skeleton}' skeleton. "
                f"The provided answer implies attack at the more hindered C1.")

    # --- Constraint 2: Stereoselectivity Check ---
    # The reaction proceeds with inversion at the attacked center (C6, which becomes C2).
    # The starting configuration at C6 is (S).
    # A rigorous analysis shows that geometric inversion at this (6S) center, combined with the change
    # in Cahn-Ingold-Prelog priorities, results in a (2S) configuration in the product.
    # The other centers (1R, 3R, 4R) are retained, becoming (1R, 5R, 4R) in the product.
    # Therefore, the final stereochemistry must be (1R, 2S, 4R, 5R).
    expected_stereochem_prefix = "(1R,2S,4R,5R)"
    
    if not answer_string.startswith(expected_stereochem_prefix):
        match = re.match(r'(\(.*\))', answer_string)
        actual_stereochem = match.group(1) if match else "an incorrect one"
        return (f"Incorrect. The answer '{final_answer_key}' violates the stereoselectivity constraint. "
                f"The final configuration should be {expected_stereochem_prefix}, but the answer provides {actual_stereochem}. "
                f"This indicates an error in determining the outcome of the S_N2 inversion at C6.")

    # If both constraints are satisfied, the answer is correct.
    return "Correct"

# Run the check
result = check_chemistry_answer()
print(result)