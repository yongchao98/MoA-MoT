import sympy

def get_prime_factors(n):
    """
    Returns the set of prime factors of an integer n.
    """
    n = abs(n)
    if n == 0 or n == 1:
        return set()
    # sympy.factorint returns a dict of prime factors and their powers
    return set(sympy.factorint(n).keys())

def analyze_curve(name, f):
    """
    Analyzes a curve z^2 = f(x) for good reduction properties above 2.
    """
    print(f"Analyzing curve {name}: z^2 = {f}")
    
    # Use sympy.Poly to handle polynomial properties like the leading coefficient
    poly_f = sympy.Poly(f)
    lc = poly_f.LC()
    
    # Calculate the discriminant of the polynomial
    disc = sympy.discriminant(f)
    
    # Find prime factors of the discriminant and the leading coefficient
    disc_factors = get_prime_factors(disc)
    lc_factors = get_prime_factors(lc)
    
    # Bad reduction for an odd prime p occurs if p divides the discriminant or the leading coefficient.
    bad_reduction_primes = {p for p in disc_factors if p != 2}.union({p for p in lc_factors if p != 2})

    print(f"  Leading Coefficient: {lc}")
    print(f"  Discriminant: {disc}")
    # Display the discriminant as a product of its prime factors
    print(f"  Factored Discriminant: {sympy.factorint(disc)}")
    
    if not bad_reduction_primes:
        print("  Result: This curve has good reduction for all odd primes (p > 2).")
    else:
        print(f"  Result: This curve has bad reduction at odd prime(s): {sorted(list(bad_reduction_primes))}")
    print("-" * 50)

# Define the symbol for our polynomials
x = sympy.symbols('x')

# Define the polynomials from each answer choice
curves = {
    'A': x**5 + 3,
    'B': x**5 - 1,
    'C': x**6 - 1,
    'D': 2*x**5 + 2*x**3 + 1,
    'E': 4*x + x**2 + 4*x**3 + 4*x**5
}

# Reorder E to standard form to be explicit, though sympy handles it
curves['E'] = 4*x**5 + 4*x**3 + x**2 + 4*x

# Analyze each curve
for name, f_poly in curves.items():
    analyze_curve(name, f_poly)