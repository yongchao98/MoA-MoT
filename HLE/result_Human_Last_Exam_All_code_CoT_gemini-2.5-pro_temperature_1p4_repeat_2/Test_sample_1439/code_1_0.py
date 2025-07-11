import sympy

def demonstrate_nu_expansion():
    """
    Demonstrates the series expansion of the critical exponent nu
    in the coupling constant u.
    """
    # Define symbols for the coupling constant 'u' and the coefficients
    # of the expansion of the anomalous dimension gamma_m.
    # From one-loop calculations, the series for gamma_m starts with a linear term in u.
    # gamma_m(u) = A*u + B*u^2 + O(u^3)
    u = sympy.Symbol('u')
    A = sympy.Symbol('A')
    B = sympy.Symbol('B')

    # The anomalous dimension gamma_m(u) as a series in u.
    # The first term is of order u, which is non-vanishing (A != 0).
    gamma_m = A * u + B * u**2

    # The critical exponent nu is a function of gamma_m.
    # nu(u) = 1 / (2 - gamma_m(u))
    nu_of_u = 1 / (2 - gamma_m)

    # Perform a Taylor series expansion of nu(u) around u=0 to see the
    # structure of the corrections. We expand up to the u^2 term.
    nu_series = nu_of_u.series(u, 0, 3)

    print("The critical exponent nu is related to the coupling constant u by the formula:")
    print("ν(u) = 1 / (2 - γ_m(u))")
    print("\nThe perturbative expansion of the anomalous dimension γ_m(u) starts at order u:")
    print(f"γ_m(u) = {gamma_m} + O(u³)")
    print("\nExpanding ν(u) as a Taylor series around u = 0 gives:")
    print(f"ν(u) = {nu_series}")

    # Extract the terms of the series to analyze the order
    zero_order_term = nu_series.coeff(u, 0)
    first_order_term_coeff = nu_series.coeff(u, 1)

    print("\nAnalysis of the expansion:")
    print(f"The 0th order term (mean-field value) is: {zero_order_term}")
    print(f"The coefficient of the u^1 term is: {first_order_term_coeff}")
    print("\nSince the coefficient 'A' is non-zero, the first non-vanishing contribution to ν")
    print("beyond the mean-field value is of the first order in the coupling constant u.")

demonstrate_nu_expansion()