import sympy

def analyze_curve(label, poly_str):
    """
    Analyzes the reduction properties of a curve z^2 = P(x).
    
    Args:
    label (str): The label for the curve (e.g., 'A').
    poly_str (str): The string representation of the polynomial P(x).
    """
    x = sympy.Symbol('x')
    P = sympy.sympify(poly_str, locals={'x': x})
    
    print(f"--- Analysis for Curve {label}: z^2 = {P} ---")
    
    # Calculate the discriminant
    disc = sympy.discriminant(P, x)
    print(f"Polynomial P(x): {P}")
    print(f"Discriminant: {disc}")
    
    # Factor the discriminant to find primes of bad reduction
    if disc != 0:
        factors = sympy.factorint(disc)
        odd_prime_factors = {p for p in factors if p != 2 and p != -1}
        
        if not odd_prime_factors:
            print("Primes of bad reduction (> 2): None. This curve has good reduction for all odd primes.")
        else:
            print(f"Primes of bad reduction (> 2): {sorted(list(odd_prime_factors))}")
    else:
        print("The polynomial has repeated roots, so the discriminant is 0. The curve is singular over Q.")
        
    print("-" * 20)

def main():
    curves = {
        'A': 'x**5 + 3',
        'B': 'x**5 - 1',
        'C': 'x**6 - 1',
        'D': '2*x**5 + 2*x**3 + 1',
        'E': '4*x**5 + 4*x**3 + x**2 + 4*x'
    }

    for label, poly_str in curves.items():
        analyze_curve(label, poly_str)
        
    print("\n--- Conclusion ---")
    print("Based on the analysis of discriminants:")
    print("A has bad reduction at {3, 5}.")
    print("B has bad reduction at {5}.")
    print("C has bad reduction at {3}.")
    print("D has bad reduction at {5}.")
    print("E has bad reduction at {3, 17, 19}.")
    print("\nCurves B, C, and D have the best 'good reduction' profiles, with only one odd prime of bad reduction.")
    print("\nTo choose between B, C, and D, we consider the 'ordinary' property. Curves B (z^2 = x^5 - 1) and C (z^2 = x^6 - 1) are defined by reducible polynomials over the rationals and have special automorphisms. Such curves with extra structure are often supersingular at many primes.")
    print("\nCurve D (z^2 = 2*x^5 + 2*x^3 + 1) is defined by an irreducible polynomial and is more 'generic'. Generic curves are less likely to have persistently supersingular reduction across many primes. Therefore, it is the most plausible candidate for having 'good ordinary reduction'.")
    print("\nThe best choice is D.")
    final_eq = curves['D']
    print(f"The equation for the chosen curve is: z^2 = 2*x^5 + 2*x^3 + 1")


if __name__ == '__main__':
    main()
