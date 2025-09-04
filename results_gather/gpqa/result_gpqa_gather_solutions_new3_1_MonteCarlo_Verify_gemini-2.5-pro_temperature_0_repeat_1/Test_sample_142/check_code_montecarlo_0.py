import re

def check_answer_correctness(final_answer_str):
    """
    Checks the correctness of the final answer for the Pinacol rearrangement question.

    This function encodes the chemical principles to determine the correct starting material 'A' 
    and product 'B', and then verifies if the chosen option matches this correct combination.

    Args:
        final_answer_str (str): The final answer string from the LLM, e.g., "<<<D>>>".

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the reason.
    """

    # --- Step 1: Define the options as provided in the original question ---
    options = {
        'A': {
            "A": "1-(hydroxydi-p-tolylmethyl)cyclohexan-1-ol",
            "B": "methyl 2-methyl-3-oxo-2-(p-tolyl)propanoate"
        },
        'B': {
            "A": "1-(hydroxydi-p-tolylmethyl)cyclopentan-1-ol",
            "B": "methyl 2-methyl-3-oxo-2-(p-tolyl)propanoate"
        },
        'C': {
            "A": "1-(hydroxydi-p-tolylmethyl)cyclohexan-1-ol",
            "B": "methyl 3-oxo-2-(p-tolyl)butanoate"
        },
        'D': {
            "A": "1-(hydroxydi-p-tolylmethyl)cyclopentan-1-ol",
            "B": "methyl 3-oxo-2-(p-tolyl)butanoate"
        }
    }

    # --- Step 2: Determine the correct answer based on chemical principles ---

    # Principle for Reaction 1: Ring expansion from a 5-membered to a 6-membered ring is favorable.
    correct_A_name = "1-(hydroxydi-p-tolylmethyl)cyclopentan-1-ol"

    # Principle for Reaction 2: Formation of the most stable carbocation (tertiary, benzylic) followed by a 1,2-hydride shift.
    correct_B_name = "methyl 3-oxo-2-(p-tolyl)butanoate"

    # Find which option letter matches the correct combination of A and B
    correct_option_letter = None
    for letter, components in options.items():
        if components["A"] == correct_A_name and components["B"] == correct_B_name:
            correct_option_letter = letter
            break
    
    if correct_option_letter is None:
        # This is a sanity check for the script itself.
        return "Error in checker: Could not find an option that matches the derived correct answer based on chemical principles."

    # --- Step 3: Check the provided final answer against the derived correct answer ---
    try:
        # Extract the letter from the final answer string
        match = re.search(r'<<<([A-D])>>>', final_answer_str)
        if not match:
            return f"Invalid answer format. Expected '<<<X>>>' where X is A, B, C, or D. Got: {final_answer_str}"
        
        chosen_letter = match.group(1)
        
        if chosen_letter == correct_option_letter:
            return "Correct"
        else:
            # Provide a detailed reason why the chosen answer is wrong.
            chosen_A = options[chosen_letter]["A"]
            chosen_B = options[chosen_letter]["B"]
            
            reason = f"The answer '{chosen_letter}' is incorrect. The correct answer is '{correct_option_letter}'.\n"
            
            if chosen_A != correct_A_name:
                reason += (
                    f"Reason for 'A': The starting material 'A' is incorrect. "
                    f"The reaction produces a six-membered ring (cyclohexanone). This is the result of a favorable ring-expansion from a five-membered ring starting material ('{correct_A_name}'). "
                    f"The chosen answer incorrectly identifies 'A' as '{chosen_A}', which would lead to a seven-membered ring product.\n"
                )
            
            if chosen_B != correct_B_name:
                reason += (
                    f"Reason for 'B': The product 'B' is incorrect. "
                    f"The reaction proceeds via the most stable carbocation (at C2) followed by a 1,2-hydride shift (H has higher migratory aptitude than methyl). "
                    f"This mechanism leads to the formation of a ketone at C3, making the product '{correct_B_name}'. "
                    f"The chosen answer incorrectly identifies 'B' as '{chosen_B}'."
                )
            
            return reason.strip()

    except Exception as e:
        return f"An error occurred while checking the answer: {e}"

# The final answer to be checked is "<<<D>>>"
# print(check_answer_correctness("<<<D>>>"))