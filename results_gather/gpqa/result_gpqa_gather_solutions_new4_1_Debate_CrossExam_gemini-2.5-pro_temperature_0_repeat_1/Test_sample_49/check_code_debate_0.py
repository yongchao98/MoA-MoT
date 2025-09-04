import re

def check_feynman_loop_answer():
    """
    Checks the correctness of the final answer for the Feynman diagram loop question.
    
    The logic is based on the standard rule in 4D Quantum Field Theory:
    - Each loop (L) contributes a factor of 1/(4pi)^2.
    - L loops contribute a factor of (1/(4pi)^2)^L = 1/(4pi)^(2L).
    
    The code will:
    1. Parse the mathematical expression to find the exponent of (4pi).
    2. Calculate the number of loops based on the rule 2L = exponent.
    3. Parse the provided final answer to get the chosen option.
    4. Compare the calculated number of loops with the value of the chosen option.
    """
    
    # --- Problem and Answer Data ---
    
    # The expression from the question
    expression = "alpha^3 * g^2 sqrt(2) * 8 * 1/(4pi)^6 * (Q/M)^2"
    
    # The options provided in the question
    options = {'A': 3, 'B': 1, 'C': 6, 'D': 2}
    
    # The final answer provided by the LLM to be checked
    final_answer = "<<<A>>>"
    
    # --- Verification Logic ---
    
    # 1. Find the exponent of (4pi) in the expression
    match = re.search(r"1/\(4pi\)\^(\d+)", expression)
    
    if not match:
        return "Failure in Check: Could not find the characteristic loop factor '1/(4pi)^n' in the expression."
        
    exponent = int(match.group(1))
    
    # 2. Calculate the number of loops (L)
    # The rule is 2 * L = exponent
    if exponent % 2 != 0:
        return f"Incorrect: The exponent in the loop factor is {exponent}, which is odd. This violates the rule 2*L = exponent, suggesting a misunderstanding of the problem's physics."
        
    calculated_loops = int(exponent / 2)
    
    # 3. Extract the chosen option from the final answer string
    answer_choice_match = re.search(r"<<<([A-D])>>>", final_answer)
    
    if not answer_choice_match:
        return f"Failure in Check: The final answer format '{final_answer}' is invalid."
        
    chosen_option_letter = answer_choice_match.group(1)
    
    # 4. Get the numerical value corresponding to the chosen option
    chosen_option_value = options.get(chosen_option_letter)
    
    # 5. Compare the calculated result with the chosen answer
    if calculated_loops == chosen_option_value:
        return "Correct"
    else:
        return (f"Incorrect: The final answer is inconsistent with the physical principles. "
                f"The expression contains the term '1/(4pi)^{exponent}', which implies 2*L = {exponent}. "
                f"Therefore, the number of loops must be L = {calculated_loops}. "
                f"The provided answer was option {chosen_option_letter}, which corresponds to {chosen_option_value} loops. "
                f"The calculated value ({calculated_loops}) does not match the answer's value ({chosen_option_value}).")

# Execute the check and print the result
result = check_feynman_loop_answer()
print(result)