import re

def check_feynman_loop_answer(proposed_answer: str):
    """
    Checks the correctness of the proposed answer for the number of loops in a Feynman diagram.

    The number of loops (L) in a Feynman diagram calculation in 4D spacetime
    is determined by the power of the loop factor 1/(4pi)^2. The total factor
    from L loops is (1/(4pi)^2)^L = 1/(4pi)^(2L).

    Args:
        proposed_answer (str): The letter choice for the answer (e.g., 'A', 'B', 'C', 'D').

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """
    # The expression given in the problem
    expression = "alpha^3 * g^2 sqrt(2) * 8 * 1/(4pi)^6 * (Q/M)^2"

    # Map the answer choices to their corresponding numerical values (number of loops)
    answer_mapping = {
        'A': 6,
        'B': 2,
        'C': 3,
        'D': 1
    }

    if proposed_answer not in answer_mapping:
        return f"Invalid answer choice '{proposed_answer}'. The choice must be one of {list(answer_mapping.keys())}."

    proposed_loops = answer_mapping[proposed_answer]

    # Use regex to find the term '1/(4pi)^n' and extract the exponent 'n'
    match = re.search(r'1/\(4pi\)\^(\d+)', expression)

    if not match:
        return "Constraint not satisfied: The expression string does not contain the expected loop factor '1/(4pi)^n'."

    # The exponent of (4pi) in the denominator
    power_of_4pi = int(match.group(1))

    # The physical rule relating the exponent to the number of loops (L) is: 2 * L = power_of_4pi
    # Therefore, L = power_of_4pi / 2
    if power_of_4pi % 2 != 0:
        return (f"Constraint not satisfied: The power of (4pi) is {power_of_4pi}, which is odd. "
                f"According to loop calculation conventions, this power should be an even number (2 * number of loops).")

    calculated_loops = power_of_4pi / 2

    # Check if the calculated number of loops matches the proposed answer's number of loops
    if calculated_loops == proposed_loops:
        return "Correct"
    else:
        reason = (
            f"Incorrect. The proposed answer '{proposed_answer}' corresponds to {proposed_loops} loops.\n"
            f"The correct number of loops is derived from the expression's loop factor.\n"
            f"Rule: Each loop in 4D spacetime contributes a factor of 1/(4pi)^2.\n"
            f"The expression contains the term '1/(4pi)^{power_of_4pi}'.\n"
            f"Therefore, the number of loops (L) is found by solving 2 * L = {power_of_4pi}.\n"
            f"This gives L = {int(calculated_loops)} loops."
        )
        return reason

# To use the checker, you would pass the letter of the answer you want to verify.
# For example, to check if 'C' is the correct answer:
# print(check_feynman_loop_answer('C'))
#
# To check if 'A' is the correct answer:
# print(check_feynman_loop_answer('A'))
