import sympy

def analyze_curve_reduction(label, poly_str):
    """
    Analyzes the reduction properties of a hyperelliptic curve z^2 = P(x)
    by calculating the prime factors of the discriminant of P(x).
    """
    x = sympy.Symbol('x')
    try:
        P = sympy.sympify(poly_str)
        disc = sympy.discriminant(P, x)
    except Exception as e:
        print(f"Could not process curve {label}: {poly_str}\nError: {e}")
        return

    print(f"Analyzing Curve {label}: z^2 = {poly_str}")

    if disc == 0:
        print("  - The discriminant is 0, so the polynomial has repeated roots.")
        print("  - The curve is singular over the rational numbers.")
        print("-" * 30)
        return

    # Factor the integer discriminant to find primes of bad reduction
    try:
        factorization = sympy.factorint(disc)
    except Exception as e:
        print(f"  - Discriminant is {disc}, but it could not be factored: {e}")
        print("-" * 30)
        return
        
    # Build the equation string showing the discriminant and its factorization
    # The final equation is Discriminant = p1^e1 * p2^e2 * ...
    equation_parts = [f"Discriminant = {disc}"]
    
    # Handle negative sign
    display_factors = []
    if -1 in factorization:
        display_factors.append("(-1)")
        del factorization[-1]
        
    display_factors.extend([f"{p}^{e}" for p, e in sorted(factorization.items())])
    equation_parts.append(" = ".join(display_factors))

    print(f"  - {' = '.join(equation_parts)}")

    # Identify odd prime factors
    bad_primes = [p for p in factorization if p > 2]

    if not bad_primes:
        print("  - This curve has good reduction for all primes p > 2.")
    else:
        print(f"  - This curve has bad reduction at prime(s): {sorted(bad_primes)}")
    
    print("-" * 30)

if __name__ == '__main__':
    curves = {
        "A": "x**5 + 3",
        "B": "x**5 - 1",
        "C": "x**6 - 1",
        "D": "2*x**5 + 2*x**3 + 1",
        "E": "4*x + x**2 + 4*x**3 + 4*x**5"
    }

    for label, poly_str in curves.items():
        analyze_curve_reduction(label, poly_str)