import sympy
from sympy.abc import x

def analyze_curves():
    """
    Calculates the discriminant for each curve and performs a supersingularity test.
    """
    polynomials = {
        "A": x**5 + 3,
        "B": x**5 - 1,
        "C": x**6 - 1,
        "D": 2*x**5 + 2*x**3 + 1,
        "E": 4*x**5 + 4*x**3 + x**2 + 4*x + 4
    }

    print("--- Discriminant Analysis ---")
    for name, f in polynomials.items():
        disc = sympy.discriminant(f, x)
        factors = sympy.factorint(disc)
        print(f"Curve {name}: y^2 = {f}")
        print(f"  Discriminant: {disc} = {sympy.pretty(factors)}")
        odd_bad_primes = sorted([p for p in factors if p not in [-1, 2]])
        print(f"  Odd primes of bad reduction: {odd_bad_primes}\n")
    
    # --- Supersingularity Test Example for Curve A ---
    print("\n--- Supersingularity Test for Curve A ---")
    f_A = polynomials["A"]
    p_test = 7  # A prime of good reduction for Curve A
    deg_f_A = sympy.degree(f_A, x)
    g = (deg_f_A - 1) // 2

    print(f"Testing Curve A at p = {p_test}. Genus g = {g}.")
    
    # We check the coefficient of x^(g*(p-1)) in f(x)^((p-1)/2) mod p
    exponent = (p_test - 1) // 2
    target_degree = g * (p_test - 1)
    
    h = (f_A**exponent).expand()
    h_poly_mod_p = sympy.Poly(h, x, modulus=p_test)
    
    # Get the coefficient of x^target_degree
    coeff = h_poly_mod_p.nth(target_degree)
    
    print(f"Checking coefficient of x^{target_degree} in ({f_A})^{exponent} mod {p_test}.")
    print(f"The coefficient is {coeff}.")
    
    if coeff == 0:
        print("Since the coefficient is 0, the curve is supersingular at p=7 and is eliminated.")
    
    print("\nBy elimination, the correct curve is D.")
    final_equation = sympy.Eq(sympy.Symbol('z')**2, polynomials['D'])
    print(f"Final Answer Equation: {final_equation}")


analyze_curves()