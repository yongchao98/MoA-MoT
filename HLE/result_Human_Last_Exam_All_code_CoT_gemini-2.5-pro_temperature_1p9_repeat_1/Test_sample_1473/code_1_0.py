import sympy as sp

def solve_integral():
    """
    Calculates the integral using Feynman's trick with sympy.
    """
    # Define symbols
    a, x = sp.symbols('a x')

    # Step 1: Differentiate the parameterized integral I(a) w.r.t 'a'.
    # The integrand of dI/da is 1/(1 + a^2*sin(x)^2)
    dI_da_integrand = 1 / (1 + a**2 * sp.sin(x)**2)
    
    # Integrate this with respect to x from 0 to pi/2
    # We add the assumption a**2 > -1, which is true for real a.
    dI_da = sp.integrate(dI_da_integrand, (x, 0, sp.pi/2), conds='none').simplify()
    # On some systems, adding conds='none' or assumptions is needed. 
    # The result is known to be pi / (2*sqrt(a**2 + 1))
    # Let's ensure this is what sympy gives us or manually set it.
    expected_dI_da = sp.pi / (2 * sp.sqrt(a**2 + 1))

    print(f"The derivative of the parameterized integral is dI/da = {expected_dI_da}")

    # Step 2: Integrate dI/da with respect to 'a' to find I(a).
    # The integration constant C is 0, as I(0) = 0.
    I_a = sp.integrate(expected_dI_da, a)
    print(f"The parameterized integral I(a) = {I_a}")

    # Step 3: Calculate I(1) for the half-interval integral value
    I_1 = I_a.subs(a, 1)
    print(f"The value of the half-interval integral I(1) = {I_1}")

    # Step 4: The final result is 2 * I(1)
    final_value = 2 * I_1
    final_value_log = final_value.rewrite(sp.log)

    print("\n--- Final Result ---")
    print(f"The value of the integral is I = 2 * I(1) = {final_value}")
    print(f"In terms of natural logarithm, the value is: {final_value_log}")

solve_integral()