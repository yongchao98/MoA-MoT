import sympy

def derive_nu_correction_order():
    """
    Symbolically derives the order of the first correction to the critical exponent nu.
    """

    # Define symbols for the derivation.
    # u_star is the coupling constant 'u' at the non-trivial fixed point.
    u_star = sympy.Symbol('u*')
    # A is a positive constant from the one-loop calculation of eta_r.
    A = sympy.Symbol('A', positive=True)

    print("Step 1: The Mean-Field (zeroth-order) value for the critical exponent nu is:")
    nu_mf = sympy.Rational(1, 2)
    print(f"ν_mf = {nu_mf}\n")

    print("Step 2: The exact relation for nu from the Renormalization Group is:")
    # eta_r is the anomalous dimension of the phi^2 operator.
    eta_r = sympy.Symbol('η_r')
    nu_expr = 1 / (2 + eta_r)
    print(f"ν = 1 / (2 + η_r)\n")

    print("Step 3: At one-loop order, eta_r is linear in the coupling constant u.")
    print("We evaluate this at the non-trivial fixed point u*:")
    # Substitute the one-loop expression for eta_r.
    eta_r_one_loop = -A * u_star
    print(f"η_r(u*) = {eta_r_one_loop} + O(u*²)\n")

    print("Step 4: Substitute the one-loop expression for eta_r into the formula for nu:")
    nu_at_u_star = 1 / (2 + eta_r_one_loop)
    print(f"ν(u*) = 1 / (2 + ({eta_r_one_loop}))")
    # Simplify the expression
    nu_at_u_star_simplified = 1 / (2 - A * u_star)
    print(f"      = {nu_at_u_star_simplified}\n")


    print("Step 5: Perform a Taylor expansion in u* to find the first correction.")
    # Series expansion around u* = 0 to see the correction terms.
    nu_series = nu_at_u_star_simplified.series(u_star, 0, 3)
    print("Expanding for small u* gives:")
    print(f"ν(u*) = {nu_series}\n")

    # Extract the terms from the series expansion
    nu_0 = nu_series.coeff(u_star, 0)
    nu_1_coeff = nu_series.coeff(u_star, 1)

    print("Conclusion:")
    print("The final equation for ν, expanded to the first order in u*, is:")
    print(f"ν = {nu_0} + ({nu_1_coeff}) * u* + O(u*²)")
    print("\nThe mean-field value is 1/2. The first non-vanishing contribution is the term")
    print(f"proportional to u*, which is of the first order in the coupling constant.")

if __name__ == '__main__':
    derive_nu_correction_order()
