import sympy
from sympy import factorint

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    root = sympy.isqrt(n)
    return root * root == n

def analyze_curve(label, poly_str):
    """
    Analyzes a curve of the form z^2 = f(x) by computing the discriminant of f(x).
    """
    x = sympy.symbols('x')
    try:
        # Create a polynomial object from the string
        f = sympy.Poly(poly_str, x)
        
        # Calculate the discriminant
        disc = sympy.discriminant(f)
        
        # Get the prime factorization of the discriminant
        # We take the absolute value for factorization
        factors = factorint(abs(disc))
        
        # Check for odd prime factors
        odd_prime_factors = {p for p in factors if p != 2}
        
        # Check if the discriminant is a perfect square
        is_sq = is_perfect_square(disc)
        
        print(f"Curve {label}: z^2 = {f.as_expr()}")
        print(f"  Discriminant: {disc}")
        print(f"  Prime Factorization of |Discriminant|: {factors}")
        
        if not odd_prime_factors:
            print("  This curve has good reduction for all odd primes.")
        else:
            print(f"  This curve has bad reduction at odd primes: {sorted(list(odd_prime_factors))}.")
            
        if is_sq:
            print(f"  A special property: The discriminant is a perfect square.")
        print("-" * 20)

    except Exception as e:
        print(f"Could not process curve {label}: {e}")

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
