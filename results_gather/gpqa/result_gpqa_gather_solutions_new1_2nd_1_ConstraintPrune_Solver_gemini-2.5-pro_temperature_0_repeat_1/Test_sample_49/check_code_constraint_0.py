import re

def check_feynman_loop_answer():
    """
    Checks the correctness of the answer to the Feynman diagram loop question.
    """
    # --- Problem Data ---
    # The expression from the colleague's note
    expression = "alpha^3 * g^2 sqrt(2) * 8 * 1/(4pi)^6 * (Q/M)^2"
    
    # The multiple-choice options provided in the question
    options = {
        'A': 6,
        'B': 2,
        'C': 1,
        'D': 3
    }
    
    # The final answer provided by the LLM
    llm_final_answer_letter = 'D'

    # --- Verification Logic ---
    
    # 1. Extract the exponent from the loop factor term '1/(4pi)^...'
    # This is the most crucial piece of information.
    match = re.search(r'1/\(4pi\)\^(\d+)', expression)
    
    if not match:
        return "Incorrect: The checker code could not find the loop factor '1/(4pi)^...' in the expression."
        
    total_exponent = int(match.group(1))

    # 2. Apply the physical principle
    # In 4D spacetime, each loop (L) contributes a factor of 1/(4pi)^2.
    # So, for L loops, the total exponent is 2 * L.
    # We need to solve: 2 * L = total_exponent
    
    # The exponent must be an even number for the principle to hold directly.
    if total_exponent % 2 != 0:
        return (f"Incorrect: The exponent in the loop factor ({total_exponent}) is not an even number, "
                f"which contradicts the standard physical principle that each loop contributes a factor of 1/(4pi)^2.")

    calculated_loops = total_exponent / 2

    # 3. Check if the LLM's answer matches the calculated result
    if llm_final_answer_letter not in options:
        return f"Incorrect: The provided answer '{llm_final_answer_letter}' is not a valid option choice."

    llm_answer_value = options[llm_final_answer_letter]

    if calculated_loops == llm_answer_value:
        return "Correct"
    else:
        # Find the correct option letter for the calculated result
        correct_letter = [key for key, val in options.items() if val == calculated_loops]
        
        reason = (
            f"Incorrect: The provided answer is wrong.\n"
            f"Reasoning: The expression contains the term '1/(4pi)^{total_exponent}'. "
            f"According to standard quantum field theory in 4D spacetime, each loop contributes a factor of '1/(4pi)^2'. "
            f"Therefore, the number of loops (L) is determined by the equation 2 * L = {total_exponent}, which gives L = {int(calculated_loops)}. "
            f"The correct number of loops is {int(calculated_loops)}, which corresponds to option {correct_letter[0]}.\n"
            f"The provided answer was option {llm_final_answer_letter}, which corresponds to {llm_answer_value} loops."
        )
        return reason

# Run the check and print the result
result = check_feynman_loop_answer()
print(result)