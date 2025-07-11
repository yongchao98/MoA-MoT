import sympy as sp

def solve_optimal_feedback():
    """
    This function symbolically derives the optimal feedback adjustment factor 'b'
    that maximizes the channel's mutual information.
    """
    # Define the symbols used in the problem
    # P: Half of the total power budget (2P)
    # alpha: Noise correlation coefficient
    # b: Feedback adjustment factor
    # P1: Power of the first transmission
    P = sp.Symbol('P', positive=True, real=True)
    alpha = sp.Symbol('alpha', real=True)
    b = sp.Symbol('b', real=True)
    P1 = sp.Symbol('P1', positive=True, real=True)

    print("Step 1: Define the signal and noise covariance matrices (K_X and K_Z).")
    
    # The signal vector is X = [X1, X2]^T = [X1, b*Z1]^T.
    # K_X = E[X * X^T]
    # E[X1^2] = P1
    # E[X2^2] = E[(b*Z1)^2] = b^2 * E[Z1^2] = b^2 * 1 = b^2
    # E[X1*X2] = E[X1 * b*Z1] = b * E[X1] * E[Z1] = 0 (due to independence and zero mean)
    K_X = sp.Matrix([
        [P1, 0],
        [0, b**2]
    ])
    
    # The noise covariance matrix K_Z is given.
    K_Z = sp.Matrix([
        [1, alpha],
        [alpha, 1]
    ])
    
    print("Signal Covariance Matrix K_X:")
    sp.pretty_print(K_X)
    print("\nNoise Covariance Matrix K_Z:")
    sp.pretty_print(K_Z)
    
    print("\nStep 2: Apply the total power constraint tr(K_X) = 2P.")
    # The total power is tr(K_X) = P1 + b^2.
    # We set it equal to the budget 2P for maximization: P1 + b^2 = 2P.
    # From this, we can express P1 in terms of P and b.
    power_constraint_eq = sp.Eq(P1 + b**2, 2*P)
    P1_expr = sp.solve(power_constraint_eq, P1)[0]
    
    print(f"From the constraint {sp.pretty(power_constraint_eq)}, we get:")
    print(f"P1 = {P1_expr}")
    
    # Substitute P1 back into K_X
    K_X_sub = K_X.subs(P1, P1_expr)
    
    print("\nStep 3: Formulate the objective function to maximize: det(K_X + K_Z).")
    # The mutual information is maximized by maximizing det(K_X + K_Z).
    K_sum = K_X_sub + K_Z
    determinant_expr = K_sum.det()
    
    print("The expression for the determinant is:")
    sp.pretty_print(determinant_expr.simplify())
    
    print("\nStep 4: Find the optimal 'b' by maximizing the determinant.")
    # To simplify, we maximize with respect to x = b^2.
    x = sp.Symbol('x', positive=True, real=True)
    objective_func = determinant_expr.subs(b**2, x)
    
    print("Let x = b^2. The objective function becomes:")
    sp.pretty_print(objective_func)
    
    # Take the derivative with respect to x and set it to zero.
    derivative = sp.diff(objective_func, x)
    print("\nThe derivative with respect to x is:")
    sp.pretty_print(derivative)
    
    # Solve for x that makes the derivative zero.
    optimal_x_sol = sp.solve(derivative, x)
    if not optimal_x_sol:
        print("Could not find a solution.")
        return
        
    optimal_x = optimal_x_sol[0]
    
    print(f"\nSetting the derivative to zero gives the optimal value for x = b^2:")
    print(f"b^2 = {optimal_x}")
    
    # The optimal b is the square root of the optimal x.
    optimal_b = sp.sqrt(optimal_x)
    
    print("\n--- FINAL RESULT ---")
    print("The optimal feedback adjustment factor 'b' is found by taking the square root.")
    final_equation = sp.Eq(b, optimal_b)
    print("The final equation is:")
    sp.pretty_print(final_equation)

if __name__ == '__main__':
    solve_optimal_feedback()
<<<b = sqrt(P)>>>