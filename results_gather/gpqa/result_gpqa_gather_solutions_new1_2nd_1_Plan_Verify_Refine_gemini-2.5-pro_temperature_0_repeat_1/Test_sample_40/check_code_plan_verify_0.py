import sympy

def check_physics_answer():
    """
    Checks the correctness of the answer to the relativistic velocity problem.
    """
    # Define symbolic variables for the index of refraction (n) and velocity (v)
    n, v = sympy.symbols('n v', real=True, positive=True)
    
    # The speed of light in a vacuum is given as c=1
    c = 1

    # The provided options from the question
    options = {
        'A': (1 + n*v) / (n + v),
        'B': (1 - n*v) / (n + v),
        'C': (1 + n*v) / (n - v),
        'D': sympy.Integer(1)
    }

    # The final answer provided by the LLM
    llm_answer_key = 'A'
    llm_answer_expr = options.get(llm_answer_key)

    if llm_answer_expr is None:
        return f"Invalid answer key '{llm_answer_key}' provided."

    # --- Check 1: Derivation from First Principles ---
    # The speed of light in the rest frame of the glass (u') is c/n
    u_prime = c / n
    
    # The relativistic velocity addition formula for co-linear motion is: u = (u' + v) / (1 + u'v/c^2)
    derived_formula = (u_prime + v) / (1 + (u_prime * v) / c**2)
    
    # Simplify the derived formula
    simplified_derived_formula = sympy.simplify(derived_formula)
    
    # Check if the LLM's answer expression matches the derived formula
    if sympy.simplify(llm_answer_expr - simplified_derived_formula) != 0:
        return (f"Incorrect. The expression for answer {llm_answer_key} is inconsistent with the "
                f"derivation from the relativistic velocity addition formula.\n"
                f"Derived formula: {simplified_derived_formula}\n"
                f"Answer {llm_answer_key}'s formula: {llm_answer_expr}")

    # --- Check 2: Verification with Physical Constraints ---
    
    # Constraint 1: If the medium is a vacuum (n=1), the speed of light must be c=1.
    expected_result_n1 = 1
    actual_result_n1 = sympy.simplify(llm_answer_expr.subs(n, 1))
    if actual_result_n1 != expected_result_n1:
        return (f"Incorrect. The expression for answer {llm_answer_key} fails the vacuum constraint (n=1).\n"
                f"Expected result: {expected_result_n1}\n"
                f"Actual result: {actual_result_n1}")

    # Constraint 2: If the medium is stationary (v=0), the speed of light must be c/n = 1/n.
    expected_result_v0 = 1/n
    actual_result_v0 = sympy.simplify(llm_answer_expr.subs(v, 0))
    if actual_result_v0 != expected_result_v0:
        return (f"Incorrect. The expression for answer {llm_answer_key} fails the stationary medium constraint (v=0).\n"
                f"Expected result: {expected_result_v0}\n"
                f"Actual result: {actual_result_v0}")

    # Optional: Check if any other option also satisfies the constraints (uniqueness check)
    for key, expr in options.items():
        if key == llm_answer_key:
            continue
        if sympy.simplify(expr.subs(n, 1)) == expected_result_n1 and \
           sympy.simplify(expr.subs(v, 0)) == expected_result_v0:
            return (f"Ambiguous. Both answer {llm_answer_key} and {key} satisfy the physical constraints. "
                    "However, the derivation confirms A is the correct choice.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_physics_answer()
print(result)