import sympy

def analyze_curve_reduction(label, poly_str):
    """
    Analyzes the reduction properties of a curve y^2 = f(x) by factoring its discriminant.
    A curve has good reduction above 2 if its discriminant has no odd prime factors.
    """
    x = sympy.Symbol('x')
    try:
        # Parse the polynomial string
        poly = sympy.sympify(poly_str)

        # The discriminant is 0 if the polynomial has repeated roots,
        # which means the curve is singular and has bad reduction everywhere.
        if sympy.discriminant(poly, x) == 0:
            disc = 0
            factors = {}
            odd_prime_factors = "Singular curve"
        else:
            # For a non-monic polynomial p(x) = a_n*x^n + ..., the full discriminant is
            # disc(p) = a_n^(2n-2) * disc(p/a_n).
            # sympy.discriminant(poly) calculates disc(p/a_n).
            n = sympy.degree(poly, x)
            lead_coeff = poly.LC(x)
            disc = lead_coeff**(2*n - 2) * sympy.discriminant(poly, x)
            
            # Find the prime factors of the discriminant
            factors = sympy.factorint(abs(disc))
            odd_prime_factors = sorted([p for p in factors.keys() if p != 2])

        print(f"Curve {label}: z^2 = {poly}")
        print(f"  Discriminant: {disc}")
        if odd_prime_factors == "Singular curve":
            print(f"  Result: The curve is singular and does not have good reduction.")
        elif not odd_prime_factors:
            print(f"  Result: Good reduction above 2. No odd prime factors in discriminant.")
        else:
            print(f"  Result: Bad reduction at odd primes {odd_prime_factors}.")
        print("-" * 50)

    except Exception as e:
        print(f"Could not analyze curve {label}: {e}")

if __name__ == '__main__':
    # Define the polynomials for each choice
    curves = {
        'A': 'x**5 + 3',
        'B': 'x**5 - 1',
        'C': 'x**6 - 1',
        'D': '2*x**5 + 2*x**3 + 1',
        'E': '4*x**5 + 4*x**3 + x**2 + 4*x'
    }

    # Analyze each curve
    for label, poly_str in curves.items():
        analyze_curve_reduction(label, poly_str)
