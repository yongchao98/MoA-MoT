import sympy

def check_physics_answer():
    """
    Verifies the solution to the relativistic velocity addition problem using symbolic math.

    The problem asks for the speed of light (u) in a moving glass (velocity v)
    with refractive index (n), as seen by a lab observer. The speed of light in
    a vacuum (c) is taken as 1.

    The derivation follows these steps:
    1. The speed of light in the glass's rest frame is u' = c/n.
    2. The relativistic velocity addition formula is u = (u' + v) / (1 + u'*v / c**2).
    3. We substitute c=1 and u'=1/n into the formula and simplify.
    4. The simplified result is compared against the expression for option C.
    """
    try:
        # 1. Define the symbolic variables.
        # n: index of refraction
        # v: velocity of the glass
        # c: speed of light in vacuum
        n, v, c = sympy.symbols('n v c', real=True, positive=True)

        # 2. Define the speed of light in the glass's rest frame (u').
        u_prime = c / n

        # 3. Define the relativistic velocity-addition formula.
        # u is the speed of light as measured in the lab frame.
        u_lab_formula = (u_prime + v) / (1 + u_prime * v / c**2)

        # 4. Substitute the problem's constraint, c = 1.
        u_lab_with_c1 = u_lab_formula.subs(c, 1)

        # 5. Simplify the expression. sympy will substitute u_prime = 1/n and
        # perform the algebraic simplification of the complex fraction.
        derived_expression = sympy.simplify(u_lab_with_c1)

        # 6. Define the expression for the final answer provided, which is 'C'.
        # Option C is (1 + n*v) / (n + v).
        answer_c_expression = (1 + n * v) / (n + v)

        # 7. Check if the derived expression is mathematically identical to the answer's expression.
        # The `equals` method is a robust way to check for symbolic equality.
        if derived_expression.equals(answer_c_expression):
            return "Correct"
        else:
            return (f"Incorrect. The provided answer 'C' corresponds to the expression {answer_c_expression}. "
                    f"However, the correct derivation from the relativistic velocity-addition formula "
                    f"yields the expression {derived_expression}.")

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check and print the result.
result = check_physics_answer()
print(result)