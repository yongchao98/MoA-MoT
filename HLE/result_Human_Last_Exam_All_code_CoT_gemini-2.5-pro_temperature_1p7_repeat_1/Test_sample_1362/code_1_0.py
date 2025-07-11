import sympy
from sympy import Poly, discriminant, factorint
from sympy.abc import x

def analyze_curve(label, poly_expr):
    """
    Analyzes a curve z^2 = f(x) by computing and factoring the discriminant of f(x).
    """
    print(f"Curve {label}: z^2 = {poly_expr}")
    
    # Create a polynomial object
    f = Poly(poly_expr, x)
    
    # Calculate the discriminant
    disc = discriminant(f)
    print(f"  Discriminant: {disc}")
    
    # Factor the discriminant to find its prime factors
    factors = factorint(disc)
    print(f"  Prime factors of the discriminant: {factors}")
    
    odd_prime_factors = [p for p in factors if p != 2 and p != -1]
    
    if not odd_prime_factors:
        print(f"  This curve has good reduction for all primes p > 2.")
    else:
        print(f"  This curve has bad reduction at primes: {odd_prime_factors}")
    print("-" * 30)

def main():
    # Define the polynomials for each option
    curves = {
        'A': x**5 + 3,
        'B': x**5 - 1,
        'C': x**6 - 1,
        'D': 2*x**5 + 2*x**3 + 1,
        'E': 4*x**5 + 4*x**3 + x**2 + 4*x
    }
    
    for label, poly_expr in curves.items():
        analyze_curve(label, poly_expr)
        
if __name__ == '__main__':
    main()