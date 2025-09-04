import re

def check_feynman_loop_answer():
    """
    Checks the correctness of the provided answer for the Feynman loop question.

    The function verifies the answer based on the principle that in 4D QFT,
    an L-loop diagram contributes a factor of 1/(4pi)^(2L) to the amplitude.
    """
    # The expression from the problem statement
    expression = "alpha^3 * g^2 sqrt(2) * 8 * 1/(4pi)^6 * (Q/M)^2"

    # The provided answer is 'B'. The options are A) 6, B) 3, C) 2, D) 1.
    # So, the value of the provided answer is 3.
    provided_answer_value = 3

    # Use a regular expression to find the term '1/(4pi)^N' and capture the exponent N.
    # The pattern looks for '1/(4pi)^' followed by one or more digits.
    # '\(' and '\)' are used to match literal parentheses.
    match = re.search(r'1/\(4pi\)\^(\d+)', expression)

    if not match:
        return "Incorrect: The expression does not contain the characteristic loop factor '1/(4pi)^N', which is required to determine the number of loops."

    # Extract the exponent from the matched group.
    try:
        exponent = int(match.group(1))
    except (ValueError, IndexError):
        return "Incorrect: Failed to parse the exponent from the loop factor in the expression."

    # The physical relationship is: 2 * L = exponent, where L is the number of loops.
    # This implies that the exponent must be an even number.
    if exponent % 2 != 0:
        return f"Incorrect: The exponent of (4pi) is {exponent}, which is odd. This would imply a non-integer number of loops ({exponent/2}), which is physically inconsistent."

    # Calculate the number of loops
    calculated_loops = int(exponent / 2)

    # Check if the calculated number of loops matches the provided answer.
    if calculated_loops == provided_answer_value:
        return "Correct"
    else:
        reason = (f"Incorrect: The expression contains the term '1/(4pi)^{exponent}'. "
                  f"Based on the loop counting rule, this corresponds to 2L = {exponent}, "
                  f"so the number of loops L must be {calculated_loops}. "
                  f"The provided answer was {provided_answer_value}.")
        return reason

# Run the check
result = check_feynman_loop_answer()
print(result)