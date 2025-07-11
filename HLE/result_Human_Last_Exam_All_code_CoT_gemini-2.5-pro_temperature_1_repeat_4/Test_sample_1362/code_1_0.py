import sympy

def analyze_curve_reduction(label, poly_str):
    """
    Analyzes the reduction properties of a curve y^2 = f(x)
    by computing the prime factors of the discriminant of f(x).
    """
    x = sympy.Symbol('x')
    f = sympy.sympify(poly_str)

    # Calculate the discriminant of the polynomial f(x)
    disc = sympy.discriminant(f, x)

    print(f"Curve {label}: z^2 = {poly_str}")
    print(f"The discriminant of the polynomial is: {disc}")

    # A curve has bad reduction at prime p if p divides the discriminant.
    # We are looking for a curve with good reduction for all primes p > 2.
    # This means the discriminant should not have any prime factors greater than 2.
    if disc == 0:
        print("The polynomial has repeated roots, so the curve is singular over Q.")
        print("It has bad reduction at all primes.")
        return

    # Find the prime factors of the absolute value of the discriminant
    prime_factors = sympy.factorint(abs(disc))
    
    bad_reduction_primes = [p for p in prime_factors.keys() if p > 2]

    if not bad_reduction_primes:
        print("This curve has good reduction for all primes p > 2.")
    else:
        print(f"This curve has bad reduction at the following primes greater than 2: {sorted(bad_reduction_primes)}")
    print("-" * 50)


if __name__ == '__main__':
    curves = {
        'A': 'x**5 + 3',
        'B': 'x**5 - 1',
        'C': 'x**6 - 1',
        'D': '2*x**5 + 2*x**3 + 1',
        'E': '4*x + x**2 + 4*x**3 + 4*x**5'
    }

    for label, poly_str in curves.items():
        analyze_curve_reduction(label, poly_str)
