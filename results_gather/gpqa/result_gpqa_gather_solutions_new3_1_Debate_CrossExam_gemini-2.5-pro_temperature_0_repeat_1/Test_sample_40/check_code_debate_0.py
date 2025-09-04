import sympy

def check_answer():
    """
    Checks the correctness of the answer to the relativistic velocity addition problem.
    """
    # 1. Define symbolic variables for the physical quantities.
    # n: index of refraction
    # v: velocity of the glass
    n, v = sympy.symbols('n v', real=True, positive=True)
    
    # c: speed of light in vacuum, given as 1
    c = 1

    # 2. Define the known velocities based on the problem statement.
    # u_prime is the speed of light in the rest frame of the glass.
    u_prime = c / n
    
    # The velocity of the glass relative to the lab is v.

    # 3. Apply the relativistic velocity-addition formula to derive the correct answer.
    # The formula for velocities in the same direction is: u = (v_frame + v_object_in_frame) / (1 + v_frame * v_object_in_frame / c^2)
    # Here, v_frame = v and v_object_in_frame = u_prime.
    correct_expression = (v + u_prime) / (1 + (v * u_prime) / c**2)
    
    # 4. Simplify the derived correct expression.
    # sympy.simplify() will perform the algebraic manipulation, like multiplying the numerator and denominator by 'n'.
    simplified_correct_expression = sympy.simplify(correct_expression)

    # 5. Define the expression from the given answer (Option C).
    # The provided answer is C, which corresponds to the expression (1 + n*v) / (n + v).
    # Note: The order of terms in the numerator (1 + n*v) vs (n*v + 1) does not matter.
    answer_c_expression = (1 + n * v) / (n + v)

    # 6. Compare the derived correct expression with the expression from the answer.
    # We can check if the difference between the two simplified expressions is zero.
    if sympy.simplify(simplified_correct_expression - answer_c_expression) == 0:
        return "Correct"
    else:
        # This part of the code would execute if the answer were incorrect.
        reason = (f"The answer is incorrect.\n"
                  f"The correct derivation using the relativistic velocity-addition formula is:\n"
                  f"u = (v + u') / (1 + v*u'/c^2)\n"
                  f"Substituting u' = c/n = 1/n and c = 1 gives:\n"
                  f"u = (v + 1/n) / (1 + v/n)\n"
                  f"Simplifying this expression yields: {simplified_correct_expression}\n"
                  f"The expression for the given answer 'C' is: {answer_c_expression}\n"
                  f"These expressions are not equivalent.")
        return reason

# Run the check and print the result.
result = check_answer()
print(result)