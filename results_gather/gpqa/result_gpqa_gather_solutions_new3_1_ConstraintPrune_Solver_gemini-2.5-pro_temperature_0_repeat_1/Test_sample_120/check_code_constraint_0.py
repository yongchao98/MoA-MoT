import re

def check_chemistry_answer():
    """
    Checks the correctness of the final answer for the given organic chemistry question.

    The function codifies the key chemical principles (constraints) of the reaction:
    1. Regioselectivity: Attack at the less hindered carbon (C6).
    2. Stereoselectivity: Inversion of configuration at the attacked carbon (C6: S -> R).
    3. Retention of configuration at spectator centers (C1, C3, C4 all remain R).
    4. Correct IUPAC naming of the resulting product.
    """
    # The final answer provided by the LLM analysis in the prompt
    final_answer_letter = "A"

    # The options provided in the question
    options = {
        "A": "(1R,2R,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol",
        "B": "(1R,2S,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol",
        "C": "(1S,4R,5S)-2,2,4,5-tetramethylcyclohexan-1-ol",
        "D": "(1R,4R,5R)-2,2,4,5-tetramethylcyclohexan-1-ol"
    }

    chosen_answer_name = options.get(final_answer_letter)

    # --- Constraint 1: Check Regioselectivity (Product Skeleton) ---
    # Attack at less-hindered C6 gives a 1,2,4,5-substituted product.
    # Attack at more-hindered C1 would give a 2,2,4,5-substituted product.
    if "2,2," in chosen_answer_name:
        return (
            "Incorrect: The answer violates the regioselectivity constraint. "
            "The product name implies attack at the more hindered C1 carbon. "
            "Organocuprates attack the less hindered C6 carbon, which should result in a "
            "'1,2,4,5-tetramethylcyclohexan-1-ol' skeleton."
        )

    # --- Define the expected stereochemistry based on constraints ---
    expected_stereochem = {
        "1": "R",  # Retained from starting material
        "2": "R",  # Inverted from (S) at C6
        "4": "R",  # Retained from starting material
        "5": "R"   # Retained from starting material
    }

    # --- Parse the stereochemistry from the chosen answer ---
    match = re.search(r'\((.*?)\)', chosen_answer_name)
    if not match:
        return f"Error: Could not parse stereochemistry from the chosen answer: {chosen_answer_name}"
    
    chosen_stereochem_str = match.group(1)
    chosen_stereochem_parts = chosen_stereochem_str.split(',')

    # --- Check Stereochemistry Constraints ---
    for part in chosen_stereochem_parts:
        part = part.strip()
        # Extract number and configuration (e.g., '1R' -> num='1', config='R')
        num = part[:-1]
        config = part[-1]
        
        if num in expected_stereochem:
            if expected_stereochem[num] != config:
                reason = ""
                if num == '2':
                    reason = ("The S_N2 attack at C6 (which becomes the new C2) causes an "
                              "inversion of configuration from (S) to (R).")
                else:
                    reason = (f"The configuration at the spectator center C{num} should be "
                              f"retained as ({num}R).")
                
                return (
                    f"Incorrect: The stereochemistry at C{num} is wrong. "
                    f"It should be ({num}{expected_stereochem[num]}) but the answer has ({num}{config}). "
                    f"Reason: {reason}"
                )

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_chemistry_answer()
print(result)