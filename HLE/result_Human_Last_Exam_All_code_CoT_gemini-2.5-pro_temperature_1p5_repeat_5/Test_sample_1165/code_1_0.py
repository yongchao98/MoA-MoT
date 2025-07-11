import sympy

def solve_ode_fluctuations():
    """
    This script symbolically calculates the scaling of the fluctuation magnitude R(ε).
    """
    # Define symbols
    # x: spatial variable
    # L: domain length
    # epsilon: small parameter
    x, L, epsilon = sympy.symbols('x L epsilon', real=True, positive=True)

    # Step 1: Define the variance of the Green's function G0(x,z) for z ~ U(0,L).
    # Based on the derivation (integral of G0^2 minus square of integral of G0),
    # Var[G0(x,z)] = x^2 * (L-x)^2 / (12 * L^2)
    var_G0 = (x**2 * (L - x)**2) / (12 * L**2)

    # Step 2: The total variance of the fluctuation u(x) is approximately
    # Var[u(x)] ≈ ε⁴ * N * Var[G₀(x, z)].
    # Substitute L = 1/ε and N = 1/ε - 1.
    N = 1/epsilon - 1
    var_u = epsilon**4 * N * var_G0.subs(L, 1/epsilon)
    var_u = sympy.simplify(var_u)
    
    # The variance expression is:
    # var_u = ε³ * (1-ε) * x² * (1/ε - x)² / 12

    # Step 3: Find the maximum of the variance with respect to x.
    # The term x^2 * (L-x)^2 is maximized at x = L/2.
    max_var_u = var_u.subs(x, L/2).subs(L, 1/epsilon)
    max_var_u = sympy.simplify(max_var_u)
    
    # Step 4: R is the square root of the maximum variance.
    R = sympy.sqrt(max_var_u)
    R_simplified = sympy.simplify(R)
    
    # Step 5: For small epsilon, find the leading order behavior.
    R_leading_order = sympy.series(R_simplified, epsilon, 0, 2).removeO()

    print("--- Analysis of Fluctuation Magnitude R(ε) ---")
    print(f"The variance of the fluctuation u(x) is given by: Var[u(x)] = {var_u}")
    print(f"The maximum variance is found to be: max_Var = {max_var_u}")
    print("\nThe estimated maximum magnitude of fluctuations, R, is the square root of the max variance.")
    print(f"R(ε) = {R_simplified}")
    print(f"\nFor small ε, the leading order term is: R(ε) ≈ {R_leading_order}")

    # To show the numbers in the final equation: sqrt(eps)/sqrt(192) = sqrt(eps)/(8*sqrt(3))
    num, den = R_leading_order.as_numer_denom()
    den_val = den.as_coeff_Mul()[0] # Gets the constant part
    sqrt_part, int_part = sympy.perfect_power(den_val, 2)
    print(f"\nThis can be written as: R(ε) ≈ sqrt(epsilon) / ({int_part} * sqrt({sqrt_part}))")
    print("Scaling: R(ε) is proportional to epsilon^(1/2).")

    print("\n--- Analysis for z_i ~ Normal(i, 0.5) ---")
    print("For the case where z_i ~ Normal(i, 0.5), the randomness of each z_i is localized.")
    print("The variance of each term, Var[G₀(x, zᵢ)], becomes O(1) instead of O(L²) = O(ε⁻²).")
    print("The total variance thus scales as N * ε⁴ * O(1) ≈ ε⁻¹ * ε⁴ = ε³.")
    print("Therefore, R(ε) = sqrt(Var) would scale as (ε³)^(1/2) = ε^(3/2).")
    print("This scaling is different from ε^(1/2).")
    
solve_ode_fluctuations()