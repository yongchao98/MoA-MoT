import re

def check_feynman_loop_answer():
    """
    Checks the correctness of the answer for the Feynman diagram loop question.

    The core principle is that in 4D Quantum Field Theory, each loop (L)
    in a Feynman diagram contributes a factor of 1/(4pi)^2 to the final expression.
    Therefore, an L-loop diagram has a factor of (1/(4pi)^2)^L = 1/(4pi)^(2L).

    The code will:
    1. Parse the given expression to find the exponent of the (4pi) term.
    2. Use the principle 2*L = exponent to calculate the number of loops L.
    3. Map the calculated L to the given multiple-choice options.
    4. Compare the result with the provided answer.
    """
    # --- Problem Definition ---
    question_expression = "alpha^3 * g^2 sqrt(2) * 8 * 1/(4pi)^6 * (Q/M)^2"
    options = {'A': 6, 'B': 3, 'C': 2, 'D': 1}
    
    # The final answer provided by the consolidated analysis
    provided_answer_letter = "B"

    # --- Step 1: Find the exponent from the expression ---
    match = re.search(r"1/\(4pi\)\^(\d+)", question_expression)
    
    if not match:
        return "Failure in checking logic: Could not find the characteristic loop factor '1/(4pi)^N' in the expression."

    try:
        exponent = int(match.group(1))
    except (ValueError, IndexError):
        return "Failure in checking logic: Could not parse the exponent from the expression."

    # --- Step 2: Calculate the number of loops (L) ---
    # The formula is 2 * L = exponent
    if exponent % 2 != 0:
        return f"Incorrect. The exponent of (4pi) is {exponent}, which is odd. The loop formula 2*L = exponent requires an even exponent, suggesting a misinterpretation of the problem or a typo in the expression."

    calculated_loops = int(exponent / 2)

    # --- Step 3: Find the correct option letter ---
    correct_option_letter = None
    for letter, loop_count in options.items():
        if loop_count == calculated_loops:
            correct_option_letter = letter
            break
    
    if correct_option_letter is None:
        return f"Incorrect. The calculated number of loops is {calculated_loops}, which does not correspond to any of the given options: {options}."

    # --- Step 4: Compare with the provided answer and return result ---
    if provided_answer_letter == correct_option_letter:
        return "Correct"
    else:
        reason = (
            f"Incorrect. The provided answer is '{provided_answer_letter}', but the correct answer is '{correct_option_letter}'.\n"
            f"The expression contains the term '1/(4pi)^{exponent}'. "
            f"Since each loop contributes a factor of '1/(4pi)^2', the number of loops (L) is found by solving 2*L = {exponent}, which gives L = {calculated_loops}. "
            f"The number of loops is 3, which corresponds to option '{correct_option_letter}'."
        )
        return reason

# Run the check and print the result
print(check_feynman_loop_answer())