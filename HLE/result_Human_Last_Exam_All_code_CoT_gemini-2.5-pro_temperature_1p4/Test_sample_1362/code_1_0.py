import sympy

def analyze_curve(label, f_str):
    """
    Computes the discriminant of the polynomial defining a curve and finds its prime factors.
    """
    x = sympy.Symbol('x')
    f = sympy.sympify(f_str)
    
    # Ensure the input is a polynomial in x
    if not f.is_polynomial(x):
        print(f"Curve {label}: The expression is not a polynomial in x.")
        return

    poly = sympy.Poly(f, x)
    
    # Calculate the discriminant
    delta = sympy.discriminant(poly)
    
    # Find prime factors of the absolute value of the discriminant
    if delta == 0:
        prime_factors = "Discriminant is 0 (polynomial has repeated roots)"
    else:
        # factorint works on positive integers
        factors_dict = sympy.factorint(abs(delta))
        prime_factors = list(factors_dict.keys())

    print(f"Curve {label}: z^2 = {f_str}")
    print(f"  Discriminant: {delta}")
    # We are interested in odd primes of bad reduction
    odd_prime_factors = [p for p in prime_factors if p != 2]
    if not odd_prime_factors:
        print("  This curve has good reduction for all odd primes.")
    else:
        print(f"  Odd prime factors of the discriminant: {odd_prime_factors}")
    print("-" * 30)

# Analyze each curve from the answer choices
analyze_curve('A', "x**5+3")
analyze_curve('B', "x**5-1")
analyze_curve('C', "x**6-1")
analyze_curve('D', "2*x**5+2*x**3+1")
analyze_curve('E', "4*x+x**2+4*x**3+4*x**5")