import re

def check_correctness(llm_answer: str) -> str:
    """
    Checks the correctness of the LLM's answer to the Feynman diagram loop question.

    The function verifies the answer based on the following logic:
    1.  In 4D Quantum Field Theory, each loop in a Feynman diagram contributes a factor of 1/(4pi)^2 to the amplitude.
    2.  The given expression has a factor of 1/(4pi)^6.
    3.  Let L be the number of loops. The total factor is (1/(4pi)^2)^L = 1/(4pi)^(2L).
    4.  Equating the exponents gives 2L = 6, which means L = 3.
    5.  The options are A) 2, B) 6, C) 1, D) 3.
    6.  The correct numerical answer (3) corresponds to option D.
    7.  The function checks if the LLM's final answer is <<<D>>>.

    Args:
        llm_answer (str): The full text of the answer provided by the LLM,
                          which should contain the final answer in the format <<<...>>>.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the reason for the error.
    """

    # --- Step 1: Define the problem's constraints and correct solution ---

    # The expression contains 1/(4pi)^6. Each loop contributes 1/(4pi)^2.
    # The equation is 2 * L = 6.
    try:
        # Extract the exponent from the expression in the question
        # "1/(4pi)^6"
        # A more robust way would be to parse the question text, but for this specific problem,
        # hardcoding the known exponent is sufficient.
        total_exponent = 6
        exponent_per_loop = 2
        correct_loop_count = total_exponent / exponent_per_loop
    except Exception as e:
        return f"Checker Error: Could not determine the correct loop count. Details: {e}"

    if correct_loop_count != 3:
        return f"Checker Error: The calculated correct loop count is {correct_loop_count}, which is unexpected."

    # Map the options from the question to their values.
    options = {
        'A': 2,
        'B': 6,
        'C': 1,
        'D': 3
    }

    # Find the correct option letter corresponding to the calculated loop count.
    correct_option_letter = None
    for letter, value in options.items():
        if value == correct_loop_count:
            correct_option_letter = letter
            break

    if correct_option_letter is None:
        return f"Checker Error: Could not find a matching option for the correct answer (L={correct_loop_count})."

    # --- Step 2: Parse the LLM's answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not match:
        return "The answer is incorrect because the final answer format is wrong. It should be '<<<X>>>' where X is one of A, B, C, or D."

    llm_option_letter = match.group(1)

    # --- Step 3: Compare the LLM's answer with the correct solution ---
    if llm_option_letter == correct_option_letter:
        return "Correct"
    else:
        llm_chosen_value = options.get(llm_option_letter, "an invalid value")
        return (f"The answer is incorrect. "
                f"The core constraint is that each loop in 4D spacetime contributes a factor of 1/(4pi)^2. "
                f"The expression has a factor of 1/(4pi)^6. "
                f"Therefore, the number of loops L must satisfy the equation 2*L = 6, which means L = 3. "
                f"The value L=3 corresponds to option {correct_option_letter}. "
                f"The provided answer selected option {llm_option_letter}, which corresponds to {llm_chosen_value} loops, violating this constraint.")