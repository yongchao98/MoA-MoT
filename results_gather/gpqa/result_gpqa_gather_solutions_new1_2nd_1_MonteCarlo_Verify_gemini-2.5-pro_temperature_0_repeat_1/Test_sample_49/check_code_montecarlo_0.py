import re

def check_feynman_loop_answer():
    """
    Checks the correctness of the answer to the Feynman diagram loop question.

    The core principle is that in 4D spacetime, each loop (L) in a Feynman
    diagram contributes a factor of 1/(4pi)^2 to the final expression.
    Therefore, for L loops, the total factor is 1/(4pi)^(2L).
    By finding the exponent in the given expression, we can solve for L.
    """
    # --- Problem Constraints and Data ---
    # The expression from the colleague's note
    expression = "alpha^3 * g^2 sqrt(2) * 8 * 1/(4pi)^6 * (Q/M)^2"
    # The multiple-choice options provided in the question
    options = {'A': 3, 'B': 1, 'C': 6, 'D': 2}
    # The final answer provided by the LLM
    llm_answer = "<<<A>>>"

    # --- Step 1: Extract the exponent from the expression ---
    # Use regex to find the exponent in the term like '1/(4pi)^exponent'
    # We remove spaces to make the regex simpler
    match = re.search(r'1/\(4pi\)\^(\d+)', expression.replace(" ", ""))

    if not match:
        return "Incorrect: The checker could not find the loop factor '1/(4pi)^exponent' in the given expression."

    exponent = int(match.group(1))

    # --- Step 2: Calculate the number of loops (L) ---
    # The formula is 2 * L = exponent.
    # The exponent must be an even number for a valid integer number of loops.
    if exponent % 2 != 0:
        return f"Incorrect: The exponent in the loop factor is {exponent}, which is odd. The formula 2*L = exponent requires an even exponent for an integer number of loops."

    calculated_loops = int(exponent / 2)

    # --- Step 3: Extract and validate the LLM's answer ---
    answer_match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not answer_match:
        return f"Incorrect: The final answer format '{llm_answer}' is invalid. It should be in the format '<<<X>>>'."

    answer_letter = answer_match.group(1)
    if answer_letter not in options:
        return f"Incorrect: The answer '{answer_letter}' is not a valid option. Valid options are {list(options.keys())}."

    answer_value = options[answer_letter]

    # --- Step 4: Compare the calculated result with the LLM's answer ---
    if calculated_loops == answer_value:
        return "Correct"
    else:
        # Find the correct option letter for the calculated value
        correct_letter = next((letter for letter, value in options.items() if value == calculated_loops), None)
        
        reason = (
            f"Incorrect: The provided answer is {answer_letter}, which corresponds to {answer_value} loops.\n"
            f"The calculation based on the expression shows a different result.\n"
            f"1. The key term in the expression is '1/(4pi)^{exponent}'.\n"
            f"2. The standard rule in 4D QFT is that L loops contribute a factor of 1/(4pi)^(2*L).\n"
            f"3. Equating the exponents gives the equation: 2*L = {exponent}.\n"
            f"4. Solving for L gives L = {exponent} / 2 = {calculated_loops}.\n"
            f"The correct number of loops is {calculated_loops}, which corresponds to option {correct_letter}, not {answer_letter}."
        )
        return reason

# Execute the check
result = check_feynman_loop_answer()
print(result)