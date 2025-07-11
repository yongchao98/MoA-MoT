import sympy as sp

def demonstrate_gelu_approximation():
    """
    Demonstrates that a combination of two GeLU neurons can approximate a square function.

    This is done by showing that the Taylor series expansion of GeLU(a*t) + GeLU(-a*t)
    around t=0 has a leading term proportional to t^2.
    """
    # Define symbols
    t = sp.Symbol('t')
    a = sp.Symbol('a', positive=True)

    # The Gaussian Error Function (erf)
    # The standard normal CDF Phi(x) can be written using erf(x)
    # Phi(x) = 0.5 * (1 + erf(x / sqrt(2)))
    # GeLU(x) = x * Phi(x)
    
    # Define the sum of two GeLU functions with opposite inputs
    # g(t) = GeLU(a*t) + GeLU(-a*t)
    # After simplification, this becomes: g(t) = a*t * erf(a*t / sqrt(2))
    
    # We will compute the Taylor series of this expression
    term1 = a * t
    term2 = sp.erf(a * t / sp.sqrt(2))
    expression = term1 * term2

    # Compute the Taylor series of the expression around t=0 up to the 6th order
    series_expansion = sp.series(expression, t, 0, 6)

    print("The target function for each dimension is g(t) = t^2.")
    print("We can approximate this using a combination of two GeLU neurons: c * [GeLU(a*t) + GeLU(-a*t)]")
    print("\nLet's analyze the Taylor series of f(t) = GeLU(a*t) + GeLU(-a*t) around t=0:")
    print(f"f(t) = {series_expansion}\n")
    
    # Extract and print the coefficients of the series
    poly_series = series_expansion.as_poly()
    if poly_series:
        coeffs = poly_series.all_coeffs()
        # The coefficients are returned in descending order of power. Let's reverse for clarity.
        coeffs.reverse()
        terms = list(sp.ordered(poly_series.monoms()))
        terms.reverse()
        
        print("This shows the following terms:")
        # The first term is O(t^6), not a real term.
        # The actual terms are t^2 and t^4 as seen from coeffs.
        # Sympy's series handling can be tricky. Let's get coeffs directly.
        t2_coeff = series_expansion.coeff(t**2)
        t4_coeff = series_expansion.coeff(t**4)

        print(f"Term 1 (proportional to t^2): ({t2_coeff}) * t**2")
        print(f"Term 2 (proportional to t^4): ({t4_coeff}) * t**4")

        print("\nThe expansion has no constant, linear (t), or cubic (t^3) terms.")
        print("The leading term is proportional to t^2.")
        print("By choosing the constant 'c' appropriately, we can match the coefficient of t^2.")
        print("\nThis construction requires 2 neurons per dimension.")
        print("For an N-dimensional input vector, this results in a total of 2*N neurons.")

demonstrate_gelu_approximation()
