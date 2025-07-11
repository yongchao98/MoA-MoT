import sympy

def solve_critical_exponent_order():
    """
    Symbolically determines the order of the first non-vanishing contribution
    to the critical exponent nu in phi^4 theory.
    """
    # Define u as the symbolic coupling constant
    u = sympy.Symbol('u')
    # The coefficient 'A' for the one-loop contribution. It's a non-zero constant.
    A = sympy.Symbol('A')

    # The anomalous dimension gamma_m(u) starts at linear order in u.
    # gamma_m(u) = A*u + O(u^2)
    gamma_m = A * u

    # The critical exponent nu is given by nu = 1 / (2 + gamma_m(u))
    nu = 1 / (2 + gamma_m)

    # Perform a Taylor series expansion of nu around u=0 to see the contributions
    # order by order. We expand up to the O(u^2) term.
    nu_series = nu.series(u, 0, 3)

    # The mean-field value is the constant term (order u^0)
    mean_field_value = nu_series.coeff(u, 0)

    # The first correction is the term with the lowest positive power of u
    first_correction_term = nu_series.coeff(u, 1) * u

    print("In φ⁴ theory, the critical exponent ν is related to the coupling constant u by:")
    print(f"ν(u) = 1 / (2 + γ_m(u))")
    print("where the anomalous dimension γ_m(u) is expanded as a power series in u.")
    print("The first non-zero term in this expansion is linear in u: γ_m(u) = A*u + O(u²)\n")

    print("To find the first correction to the mean-field value (ν=1/2), we expand ν(u) around u=0:")
    # We use sympy.pretty_print for a nicer equation format
    print("ν(u) ≈ ", end="")
    sympy.pretty_print(nu_series)
    print("\n")

    print(f"The mean-field value (the u^0 term) is: {mean_field_value}")
    print(f"The first non-vanishing contribution (the u^1 term) is: {first_correction_term}")

    print("\nTherefore, the initial non-vanishing contribution to ν appears at the first order (u^1) in the coupling constant.")

solve_critical_exponent_order()
<<<1>>>