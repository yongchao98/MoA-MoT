import sympy as sp

def solve_critical_exponent_order():
    """
    This script symbolically derives the order in the coupling constant 'u' at which the
    critical exponent nu first receives a non-vanishing contribution in phi^4 theory.
    """
    # Define symbols for the calculation
    # u: the phi^4 coupling constant
    # epsilon: the deviation from 4 spatial dimensions (d = 4 - epsilon)
    # A, B: positive constants derived from loop integrals in the theory
    u, epsilon = sp.symbols('u epsilon')
    A, B = sp.symbols('A B', positive=True)

    print("Step-by-step derivation for the order of the first contribution to the critical exponent nu:")
    print("-" * 80)

    # --- Step 1: Define the RG functions ---
    # In the epsilon-expansion, the RG flow of the coupling 'u' is governed by the beta function.
    # To one-loop order, it is: beta(u) = -epsilon*u + B*u^2
    # The critical exponent nu is given by nu = 1 / (2 - gamma_r(u*)),
    # where gamma_r is an anomalous dimension.
    # The first contribution to gamma_r is linear in u: gamma_r(u) = A*u
    gamma_r = A * u
    
    print(f"1. The critical exponent nu is related to the coupling 'u' via an anomalous dimension gamma_r.")
    print(f"   To lowest order, gamma_r(u) is proportional to u: gamma_r(u) = {gamma_r}")

    # --- Step 2: Define nu as a function of u ---
    # The expression for nu before evaluating at the fixed point is:
    nu_of_u = 1 / (2 - gamma_r)
    
    print(f"\n2. The general expression for nu as a function of u is: nu(u) = 1 / (2 - gamma_r(u))")
    # Using sympy's pretty print for a clearer equation
    print("   Which gives:")
    sp.pprint(sp.Eq(sp.Symbol('nu(u)'), nu_of_u), use_unicode=False)
    
    # --- Step 3: Expand nu(u) as a series in u ---
    # To find the first contribution, we expand nu(u) in a power series around u=0.
    # The constant term will be the mean-field value, and the next term is the first correction.
    nu_series_in_u = nu_of_u.series(u, 0, 3)

    print("\n3. We expand nu(u) as a power series in 'u' to find the first correction:")
    # Pretty print the series expansion
    sp.pprint(sp.Eq(sp.Symbol('nu(u)'), nu_series_in_u), use_unicode=False)
    
    # --- Step 4: Identify the order ---
    print("\n4. Analysis of the series expansion:")
    print("   - The constant term is 1/2. This is the classical or mean-field value of nu.")
    print(f"   - The next term in the series is ({sp.latex(nu_series_in_u.coeff(u))})*u.")
    print("   This term is proportional to u^1.")
    print("\nConclusion:")
    print("The initial non-vanishing contribution to the critical exponent nu occurs at the first order (1) in the coupling constant u.")
    
# Execute the function
if __name__ == "__main__":
    solve_critical_exponent_order()
