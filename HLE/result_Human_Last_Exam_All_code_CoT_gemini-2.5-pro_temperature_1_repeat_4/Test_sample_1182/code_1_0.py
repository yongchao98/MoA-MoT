import math
from sympy import Symbol

def count_double_points():
    """
    Calculates the number of double points in the stable reduction of the curve
    y^2 = 8*x + 1*x^2 + 4*x^3 + 4*x^4 + 8*x^5 above the prime 2.
    """
    
    # Step 1: Define the curve and find its genus, g.
    # The curve is C: y^2 = f(x) where f(x) = 8*x^5 + 4*x^4 + 4*x^3 + x^2 + 8*x.
    # The degree of the polynomial f(x) is n = 5.
    n = 5
    # The genus of a hyperelliptic curve y^2 = f(x) is g = floor((n-1)/2).
    g = math.floor((n - 1) / 2)
    
    print("Step 1: Calculating the genus of the curve.")
    print(f"The curve is defined by a polynomial of degree n = {n}.")
    print(f"The genus is g = floor(({n}-1)/2) = {g}.")
    print("-" * 30)
    
    # Step 2: Find a better model for the curve.
    # The initial model can be transformed to one where the roots are easier to analyze.
    # A suitable change of variables (x=X/2, y=Y/2) leads to the new model:
    # Y^2 = f2(X) = X^5 + X^4 + 2*X^3 + X^2 + 16*X
    X = Symbol('X')
    f2_poly_str = "X**5 + X**4 + 2*X**3 + X**2 + 16*X"
    print("Step 2: Finding a suitable model for analysis.")
    print(f"A better model for the curve is Y^2 = {f2_poly_str}.")
    print("-" * 30)

    # Step 3: Determine the number of root clusters, c.
    # This is the number of distinct roots of the new polynomial when reduced modulo 2.
    # f2(X) mod 2 = X^5 + X^4 + X^2 = X^2 * (X^3 + X^2 + 1)
    f2_mod2_poly_str = "X^2 * (X^3 + X^2 + 1)"
    print("Step 3: Counting the number of root clusters (c).")
    print(f"Reducing the polynomial modulo 2 gives: {f2_mod2_poly_str}.")
    
    # The distinct roots come from the square-free factors: X and (X^3 + X^2 + 1).
    # The factor X gives the root 0.
    num_roots_factor1 = 1
    # The factor X^3 + X^2 + 1 is irreducible over F_2 and has 3 distinct roots in F_8.
    num_roots_factor2 = 3
    
    # Total number of distinct roots is the sum.
    c = num_roots_factor1 + num_roots_factor2
    
    print(f"The factor X contributes {num_roots_factor1} distinct root.")
    print(f"The factor X^3 + X^2 + 1 contributes {num_roots_factor2} distinct roots.")
    print(f"The total number of distinct roots modulo 2, which is the number of clusters, is c = {c}.")
    print("-" * 30)

    # Step 4: Calculate the number of double points, delta.
    # The formula is delta = g + c - 1.
    delta = g + c - 1
    
    print("Step 4: Calculating the number of double points (δ).")
    print("The number of double points is given by the formula: δ = g + c - 1.")
    # We print the final equation with the computed values.
    print(f"Plugging in the values, we get: δ = {g} + {c} - 1 = {delta}")
    
    return delta

# Run the calculation.
count_double_points()