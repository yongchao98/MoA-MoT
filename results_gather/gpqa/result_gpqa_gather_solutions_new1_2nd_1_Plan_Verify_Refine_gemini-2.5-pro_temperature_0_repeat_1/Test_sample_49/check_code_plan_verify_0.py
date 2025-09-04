import re

def check_feynman_loop_answer():
    """
    Checks the correctness of the final answer for the Feynman diagram loop question.
    The core logic is based on the principle that in 4D spacetime, L loops
    in a Feynman diagram contribute a factor of 1/(4pi)^(2L) to the final expression.
    """
    # --- Problem Constraints and Given Information ---
    # The expression from the question note
    question_expression = "alpha^3 * g^2 sqrt(2) * 8 * 1/(4pi)^6 * (Q/M)^2"
    
    # The options as stated in the final analysis. It's crucial to use the correct mapping.
    # The final answer states the options are: A) 1, B) 3, C) 6, D) 2
    options_map = {'A': 1, 'B': 3, 'C': 6, 'D': 2}
    
    # The final answer provided by the LLM to be checked.
    llm_answer = "<<<B>>>"

    # --- Step 1: Extract the exponent from the loop factor in the expression ---
    # We are looking for the 'N' in '1/(4pi)^N'.
    match = re.search(r'1/\(4pi\)\^(\d+)', question_expression)
    
    if not match:
        return "Incorrect: The checker could not find the loop factor '1/(4pi)^N' in the problem's expression."
        
    exponent = int(match.group(1))
    
    # --- Step 2: Calculate the correct number of loops (L) ---
    # The physical principle is: exponent = 2 * L
    if exponent % 2 != 0:
        return f"Incorrect: The exponent in the loop factor is {exponent}, which is odd. It must be an even number (2*L) according to QFT principles."
        
    correct_num_loops = int(exponent / 2)

    # --- Step 3: Parse the LLM's final answer ---
    try:
        # Extract the letter from the <<<...>>> format
        provided_option_letter = llm_answer.strip().replace('<', '').replace('>', '')
        if provided_option_letter not in options_map:
            return f"Incorrect: The provided answer option '{provided_option_letter}' is not a valid choice. Valid choices are {list(options_map.keys())}."
    except Exception:
        return f"Incorrect: Could not parse the final answer format: {llm_answer}"

    # Get the numerical value corresponding to the provided answer option
    provided_num_loops = options_map[provided_option_letter]

    # --- Step 4: Compare the calculated result with the provided answer ---
    if correct_num_loops == provided_num_loops:
        return "Correct"
    else:
        # Find the correct option letter for the calculated number of loops for a more descriptive error
        correct_option_letter = [key for key, val in options_map.items() if val == correct_num_loops]
        
        reason = (f"Incorrect: The calculation shows the number of loops should be {correct_num_loops}. "
                  f"The expression contains the term '1/(4pi)^6'. Based on the principle that L loops contribute a factor of '1/(4pi)^(2L)', "
                  f"we solve 2L = 6, which gives L = {correct_num_loops}. "
                  f"This corresponds to option {correct_option_letter[0]}. "
                  f"The provided answer was option {provided_option_letter}, which corresponds to {provided_num_loops} loops.")
        return reason

# Execute the check
result = check_feynman_loop_answer()
print(result)