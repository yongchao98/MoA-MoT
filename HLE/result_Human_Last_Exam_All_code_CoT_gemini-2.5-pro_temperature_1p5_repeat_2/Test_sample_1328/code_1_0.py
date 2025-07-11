import sympy

def solve_optimal_feedback_factor():
    """
    This function symbolically derives the optimal feedback adjustment factor 'b'
    that maximizes the mutual information of the described wireless channel.
    """
    # Step 1: Define symbolic variables for the problem parameters.
    # P1: Power of the first transmission
    # P2: Power of the second transmission
    # alpha: Noise correlation coefficient
    # b: Feedback adjustment factor
    P1, P2, alpha, b = sympy.symbols('P1 P2 alpha b', real=True)

    print("Plan: Maximize the determinant of the received signal's covariance matrix, |K_{X+Z}|.")
    print("-" * 70)

    # Step 2: Construct the covariance matrix of the received signal Y = X + Z.
    # The signal vector is X = [X1, X2] and the noise vector is Z = [Z1, Z2].
    # Due to feedback, the second transmission is X2 = X'_2 - b*Z1, where X'_2 is
    # the new information signal, independent of Z1.
    # The noise covariance matrix K_Z has Var(Z1)=1, Var(Z2)=1, Cov(Z1,Z2)=alpha.
    #
    # The components of the received signal Y are:
    # Y1 = X1 + Z1
    # Y2 = X2 + Z2 = (X'_2 - b*Z1) + Z2
    #
    # We compute the elements of the covariance matrix K_Y = K_{X+Z}.
    # E[X1^2] = P1, E[Z1^2] = 1, E[Z2^2] = 1, E[Z1*Z2] = alpha.
    # The power of the second transmission is specified as E[X2^2] = P2.
    # E[X2^2] = E[(X'_2 - b*Z1)^2] = E[X'_2^2] + b^2*E[Z1^2] = P'_2 + b^2.
    # Thus, the power of the new signal is P'_2 = P2 - b^2.

    # K_Y_11 = E[Y1^2] = E[(X1 + Z1)^2] = E[X1^2] + E[Z1^2] = P1 + 1
    K_Y_11 = P1 + 1

    # K_Y_22 = E[Y2^2] = E[(X'_2 - b*Z1 + Z2)^2]
    #        = E[X'_2^2] + b^2*E[Z1^2] + E[Z2^2] - 2*b*E[Z1*Z2]
    #        = (P2 - b**2) + b**2*1 + 1 - 2*b*alpha
    #        = P2 + 1 - 2*b*alpha
    K_Y_22 = P2 + 1 - 2 * b * alpha

    # K_Y_12 = E[Y1*Y2] = E[(X1 + Z1)*(X'_2 - b*Z1 + Z2)]
    #        = E[-b*Z1^2 + Z1*Z2] (other terms are zero due to independence)
    #        = -b*E[Z1^2] + E[Z1*Z2] = -b + alpha
    K_Y_12 = alpha - b

    K_Y = sympy.Matrix([[K_Y_11, K_Y_12], [K_Y_12, K_Y_22]])

    print("Step 1: The covariance matrix K_{X+Z} is constructed as:")
    print("K_{X+Z} =")
    sympy.pprint(K_Y, use_unicode=True)
    print("-" * 70)

    # Step 3: Compute the determinant of K_Y.
    det_K_Y = K_Y.det()
    det_K_Y_simplified = sympy.simplify(det_K_Y)

    print("Step 2: The determinant |K_{X+Z}| is computed as:")
    print(f"|K_{X+Z}| = {det_K_Y_simplified}")
    print("-" * 70)

    # Step 4: Maximize the determinant by finding where its derivative w.r.t. 'b' is zero.
    derivative_b = sympy.diff(det_K_Y_simplified, b)
    
    print("Step 3: To maximize, we take the derivative with respect to 'b':")
    print(f"d/db |K_{X+Z}| = {derivative_b}")
    print("-" * 70)
    
    # Solve for b to find the optimal value.
    optimal_b_solution = sympy.solve(derivative_b, b)
    optimal_b = optimal_b_solution[0]

    # Display the final result.
    final_equation = sympy.Eq(sympy.Symbol('b_opt'), optimal_b)
    print("Step 4: Solving for 'b' when the derivative is zero gives the optimal value:")
    print("Final Equation for the Optimal Feedback Factor:")
    sympy.pprint(final_equation, use_unicode=True)

if __name__ == '__main__':
    solve_optimal_feedback_factor()