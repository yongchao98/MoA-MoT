import sympy

def solve_optimal_steps():
    """
    This script calculates the optimal two-step gradient descent learning rates
    (gamma_1, gamma_2) for M-smooth and 1-strongly convex functions.

    The problem is to minimize the worst-case convergence ratio:
    ||x_2 - x_*||^2 / ||x_0 - x_*||^2

    This is equivalent to finding a polynomial P(lambda) = (1 - gamma_1*lambda)(1 - gamma_2*lambda)
    such that its maximum absolute value on the interval [1, M] is minimized.
    The solution is given by a shifted and scaled Chebyshev polynomial of degree 2.

    The optimal gamma values are the roots of the quadratic equation x^2 - S*x + P = 0,
    where S and P are the sum and product of the gammas, derived from the properties
    of the optimal polynomial.
    """

    # M is the condition number kappa
    M = sympy.Symbol('M')
    sqrt2 = sympy.sqrt(2)

    # Derived sum (S) and product (P) of the optimal step sizes (gamma_1, gamma_2)
    # S = gamma_1 + gamma_2
    S_expr = (8 * (M + 1)) / (M**2 + 6*M + 1)
    
    # P = gamma_1 * gamma_2
    P_expr = 8 / (M**2 + 6*M + 1)

    # The step sizes are the roots of x^2 - S*x + P = 0
    # Using the quadratic formula: x = (-b +- sqrt(b^2-4ac))/(2a)
    # Here, a=1, b=-S, c=P. So, gamma = (S +- sqrt(S^2-4P))/2
    discriminant = S_expr**2 - 4 * P_expr
    
    # Simplify the discriminant
    # D = (64*(M+1)**2 - 32*(M**2+6*M+1)) / (M**2+6*M+1)**2
    # D = (32*M**2 - 64*M + 32) / (M**2+6*M+1)**2
    # D = 32*(M-1)**2 / (M**2+6*M+1)**2
    # sqrt(D) = 4*sqrt(2)*(M-1) / (M**2+6*M+1)
    
    sqrt_discriminant_expr = (4 * sqrt2 * (M - 1)) / (M**2 + 6*M + 1)
    
    gamma1_expr = (S_expr - sqrt_discriminant_expr) / 2
    gamma2_expr = (S_expr + sqrt_discriminant_expr) / 2
    
    # Let's simplify the final expressions
    gamma1_expr_simplified = sympy.simplify(gamma1_expr)
    gamma2_expr_simplified = sympy.simplify(gamma2_expr)
    
    print("The problem is to find the pair of step sizes (gamma_1, gamma_2) that minimizes")
    print("the convergence rate after two steps of gradient descent.")
    print("The optimal values are derived from the roots of the Chebyshev polynomial of degree 2, scaled and shifted to the interval [1, M].\n")
    
    print("The best choice for the pair (gamma_1, gamma_2) is given by the formulas:")
    print("gamma_1 =", gamma1_expr_simplified)
    print("gamma_2 =", gamma2_expr_simplified)
    print("\nNote: The order of gamma_1 and gamma_2 can be swapped.\n")
    
    # Example calculation for a specific value of M
    M_val = 10
    print(f"--- Example Calculation for M = {M_val} ---")
    
    # Substitute M_val into expressions and show the calculation steps
    gamma1_num = gamma1_expr_simplified.subs(M, M_val)
    gamma2_num = gamma2_expr_simplified.subs(M, M_val)

    numerator1 = 4 * (M_val + 1) - 2 * sympy.sqrt(2) * (M_val - 1)
    numerator2 = 4 * (M_val + 1) + 2 * sympy.sqrt(2) * (M_val - 1)
    denominator = M_val**2 + 6 * M_val + 1
    
    print("\nStep-by-step substitution:")
    print(f"Denominator = M^2 + 6*M + 1 = {M_val}^2 + 6*{M_val} + 1 = {denominator}")
    
    print(f"\nNumerator for gamma_1 = 4*(M+1) - 2*sqrt(2)*(M-1) = 4*({M_val+1}) - 2*sqrt(2)*({M_val-1}) = {4*(M_val+1)} - {2*(M_val-1)}*sqrt(2)")
    print(f"gamma_1 = ({4*(M_val+1)} - {2*(M_val-1)}*sqrt(2)) / {denominator}")
    print(f"gamma_1 ≈ {gamma1_num.evalf()}")
    
    print(f"\nNumerator for gamma_2 = 4*(M+1) + 2*sqrt(2)*(M-1) = 4*({M_val+1}) + 2*sqrt(2)*({M_val-1}) = {4*(M_val+1)} + {2*(M_val-1)}*sqrt(2)")
    print(f"gamma_2 = ({4*(M_val+1)} + {2*(M_val-1)}*sqrt(2)) / {denominator}")
    print(f"gamma_2 ≈ {gamma2_num.evalf()}")

if __name__ == '__main__':
    solve_optimal_steps()