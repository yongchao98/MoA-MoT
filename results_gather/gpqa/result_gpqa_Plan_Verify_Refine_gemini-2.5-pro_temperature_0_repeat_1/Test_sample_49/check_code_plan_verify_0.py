import re

def check_feynman_loop_answer():
    """
    This function checks the correctness of the provided answer for the number of loops
    in a Feynman diagram based on the given mathematical expression.

    The core principle is that in 4-dimensional spacetime, an L-loop diagram's
    amplitude is generally proportional to 1/(4pi)^(2L).
    """

    # The expression from the colleague's note
    expression = "alpha^3 * g^2 sqrt(2) * 8 * 1/(4pi)^6 * (Q/M)^2"

    # The proposed answer from the LLM is C, which corresponds to 3 loops.
    llm_answer_value = 3

    # Use a regular expression to find the exponent of the (4pi) term.
    # This makes the check robust. The pattern looks for '1/(4pi)^' followed by one or more digits.
    match = re.search(r'1/\(4pi\)\^(\d+)', expression)

    # Constraint Check 1: The expression must contain the loop factor.
    if not match:
        return "Constraint not satisfied: The provided expression does not contain the expected loop factor of the form '1/(4pi)^exponent'."

    # Extract the exponent 'p' from the expression.
    p = int(match.group(1))

    # The relationship between the exponent 'p' and the number of loops 'L' is p = 2 * L.
    # We can solve for L: L = p / 2.
    calculated_loops = p / 2

    # Constraint Check 2: The number of loops must be a whole number.
    if calculated_loops != int(calculated_loops):
        return f"Constraint not satisfied: The calculation resulted in a non-integer number of loops ({calculated_loops}), which is physically impossible."

    calculated_loops = int(calculated_loops)

    # Final Check: Compare the calculated number of loops with the LLM's answer.
    if calculated_loops == llm_answer_value:
        return "Correct"
    else:
        reason = (
            f"The answer is incorrect.\n"
            f"The number of loops (L) is determined by the power of the '1/(4pi)' term.\n"
            f"The standard formula for an L-loop diagram in 4D spacetime gives a factor of 1/(4pi)^(2L).\n"
            f"The given expression has a factor of 1/(4pi)^{p}, where the exponent p is {p}.\n"
            f"By equating the exponents, we get the equation: 2L = {p}.\n"
            f"Solving for L gives L = {p} / 2 = {calculated_loops}.\n"
            f"The provided answer was {llm_answer_value}, but the correct calculation yields {calculated_loops} loops."
        )
        return reason

# Run the check
result = check_feynman_loop_answer()
print(result)