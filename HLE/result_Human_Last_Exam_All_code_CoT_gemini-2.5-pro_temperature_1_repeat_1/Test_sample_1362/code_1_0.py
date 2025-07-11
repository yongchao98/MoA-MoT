import sympy

def get_poly_and_string(curve_id):
    """Returns the polynomial f(x) and its string representation for a given curve ID."""
    x = sympy.Symbol('x')
    polys = {
        'A': x**5 + 3,
        'B': x**5 - 1,
        'C': x**6 - 1,
        'D': 2*x**5 + 2*x**3 + 1,
        'E': 4*x**5 + 4*x**3 + x**2 + 4*x
    }
    f = polys.get(curve_id)
    # Using str() for a compact representation in the printout
    curve_str = f"z^2 = {str(f)}"
    return f, curve_str

def analyze_curve(curve_id):
    """Analyzes a curve for good ordinary reduction."""
    x = sympy.Symbol('x')
    f, curve_str = get_poly_and_string(curve_id)
    print(f"--- Analyzing Curve {curve_id}: {curve_str} ---")

    # 1. Good Reduction Check (via Discriminant)
    disc = sympy.discriminant(f, x)
    try:
        factors = sympy.factorint(disc)
        bad_primes = sorted([p for p in factors if p not in [-1, 2]])
        if not bad_primes:
            print(f"Discriminant is {disc}. It has no odd prime factors.")
            print("This means the curve has good reduction for all primes p > 2.")
        else:
            print(f"Discriminant's odd prime factors are {bad_primes}.")
            print(f"This means the curve has bad reduction at p = {', '.join(map(str, bad_primes))}.")
    except Exception:
        print(f"Could not factor discriminant {disc}.")
        bad_primes = []

    # 2. Ordinary Reduction Check
    # We check at a small odd prime 'p' which is NOT a prime of bad reduction.
    primes_to_check = [3, 5, 7]
    found_failure = False
    for p in primes_to_check:
        if p in bad_primes:
            continue  # Can't check ordinarity at a prime of bad reduction

        print(f"Checking for ordinary reduction at p = {p} (a prime of good reduction)...")
        
        # Honda's criterion: check coeff of x^(p-1) in f(x)^((p-1)/2)
        k = (p - 1) // 2
        power_f = sympy.expand(f**k)
        poly_power_f = sympy.Poly(power_f, x)
        coeff = poly_power_f.nth(p - 1)
        
        print(f"The test requires the coefficient of x^{p-1} in (f(x))^{k} (where k={k}) not be zero modulo {p}.")
        print(f"The coefficient is {coeff}. The equation is {coeff} mod {p} = {coeff % p}.")

        if coeff % p == 0:
            print(f"Result: The reduction at p = {p} is NOT ordinary (it is supersingular).")
            print(f"Conclusion: Curve {curve_id} fails the test.\n")
            found_failure = True
            break # One failure is enough to disqualify the curve
        else:
            print(f"Result: The reduction at p = {p} is ORDINARY.")
            # This curve is a candidate, we break to check the next curve.
            break

    if not found_failure:
        print(f"Conclusion: Curve {curve_id} is a candidate for having good ordinary reduction above 2.\n")


# Run analysis for all curves
for curve_id in ['A', 'B', 'C', 'D', 'E']:
    analyze_curve(curve_id)
