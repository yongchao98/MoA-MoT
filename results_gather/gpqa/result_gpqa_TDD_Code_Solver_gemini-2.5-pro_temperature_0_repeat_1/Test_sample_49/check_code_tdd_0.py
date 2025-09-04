import re

def check_answer():
    """
    This function checks the correctness of the provided answer to the Feynman loop question.
    
    It verifies the logic that an L-loop diagram in 4D spacetime using dimensional
    regularization contributes a factor of 1/(4pi)^(2L) to the final expression.
    """
    
    # 1. Define the problem's input and the LLM's answer
    expression_from_question = "alpha^3 * g^2 sqrt(2) * 8 * 1/(4pi)^6 * (Q/M)^2"
    
    # The LLM's answer is C, which corresponds to 3 loops.
    llm_answer_option = "C"
    expected_loops = 3

    # 2. Apply the physical principle to derive the correct number of loops
    # We need to find the exponent 'N' in the term '1/(4pi)^N'.
    # The regex looks for the pattern '1/(4pi)^' followed by one or more digits.
    # It allows for optional whitespace around the '*' in '4*pi'.
    pattern = r'1/\(4\s*\*?\s*pi\)\^(\d+)'
    
    match = re.search(pattern, expression_from_question)
    
    if not match:
        # This case handles tree-level diagrams or expressions in a different format.
        calculated_loops = 0
    else:
        # The first captured group is the exponent.
        exponent = int(match.group(1))
        
        # The core relationship is: 2 * L = exponent
        # The exponent must be an even number for this rule to apply directly.
        if exponent % 2 != 0:
            return (f"Constraint failed: The exponent of (4pi) is {exponent}, which is an odd number. "
                    f"The standard formula L = exponent / 2 requires an even exponent.")
        
        calculated_loops = exponent // 2

    # 3. Compare the derived result with the LLM's answer
    if calculated_loops == expected_loops:
        return "Correct"
    else:
        return (f"Incorrect. The expression contains the term '1/(4pi)^6'. "
                f"According to the rule that each loop contributes a factor of '1/(4pi)^2', "
                f"the number of loops (L) is determined by solving 2L = 6. "
                f"This gives a calculated result of L = {calculated_loops} loops. "
                f"The provided answer corresponds to {expected_loops} loops, which does not match the calculation.")

# Run the check
result = check_answer()
print(result)