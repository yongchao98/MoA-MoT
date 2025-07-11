import sympy

def analyze_curve(label, poly_str):
    """
    Calculates and prints the discriminant and its prime factorization for a given polynomial.
    """
    x = sympy.Symbol('x')
    # The sympify function can parse a string into a sympy expression
    poly = sympy.sympify(poly_str, locals={'x': x})
    
    # Calculate the discriminant
    disc = sympy.discriminant(poly, x)
    
    # Factor the discriminant
    factors = sympy.factorint(disc)
    
    print(f"Curve {label}: z^2 = {poly_str}")
    print(f"  Discriminant: {disc}")
    print(f"  Prime factors of discriminant: {factors}")
    
    # Check for odd prime factors
    has_odd_prime_factor = False
    for p in factors:
        if p != 2 and p != -1:
            has_odd_prime_factor = True
            break
            
    if has_odd_prime_factor:
        print("  Result: Has bad reduction at an odd prime.\n")
    else:
        print("  Result: Has good reduction for all odd primes (p > 2).\n")


def main():
    # Polynomials from the answer choices
    curves = {
        "A": "x**5 + 3",
        "B": "x**5 - 1",
        "C": "x**6 - 1",
        "D": "2*x**5 + 2*x**3 + 1",
        "E": "4*x**5 + 4*x**3 + x**2 + 4*x"
    }

    for label, poly_str in curves.items():
        analyze_curve(label, poly_str)

if __name__ == "__main__":
    main()