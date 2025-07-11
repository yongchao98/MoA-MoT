import sympy

def solve_schwartz_moment_problem():
    """
    Explains why a Schwartz function with all zero moments must be the zero function.
    """
    print("Yes, if a Schwartz function f(x) has all its moments equal to zero, it must be the zero function.")
    print("Here is a step-by-step explanation using the Fourier transform, illustrated with sympy.\n")

    # Step 1: The Fourier Transform and its derivatives
    print("--- Step 1: The Fourier Transform and its Derivatives ---")
    print("Let F(xi) be the Fourier transform of f(x). The definition is:")
    print("F(xi) = Integral(f(x) * exp(-2*pi*i*x*xi), (x, -oo, oo))\n")
    print("The derivatives of F(xi) with respect to xi, when evaluated at xi=0, are related to the moments of f(x).")
    print("A moment is defined as M_n = Integral(x^n * f(x), (x, -oo, oo)).")
    print("Let's demonstrate the relationship for the first few derivatives (n=0, 1, 2):\n")

    x, xi = sympy.symbols('x xi')
    i, pi = sympy.I, sympy.pi
    f = sympy.Function('f')(x)

    # Loop to demonstrate the relationship for n=0, 1, 2
    for n in range(3):
        # This is the integrand of the Fourier transform definition
        integrand = f * sympy.exp(-2 * pi * i * x * xi)

        # Differentiate the integrand n times with respect to xi
        integrand_deriv_n = sympy.diff(integrand, xi, n)
        
        # Evaluate the result at xi = 0
        integrand_deriv_n_at_0 = integrand_deriv_n.subs(xi, 0)
        
        # The n-th derivative of F(xi) at xi=0 is the integral of the above
        F_deriv_n_at_0 = sympy.Integral(integrand_deriv_n_at_0, (x, -sympy.oo, sympy.oo))
        
        print(f"For n={n}:")
        print(f"The {n}-th derivative of F(xi) at xi=0, denoted F^({n})(0), is:")
        # Use pretty print for better formatting
        sympy.pprint(F_deriv_n_at_0, use_unicode=True)
        
        # Explicitly show the constant factor
        factor = (-2 * pi * i)**n
        print(f"\nThis can be rewritten as: F^({n})(0) = ({sympy.pretty(factor)}) * M_{n}")
        print("-" * 20)

    # Step 2: Applying the condition
    print("\n--- Step 2: Applying the Given Condition ---")
    print("We are given that all moments M_k are zero for k = 0, 1, 2, ...")
    print("From the relationship shown in Step 1, this directly implies that F^(k)(0) = 0 for all k >= 0.\n")

    # Step 3: Taylor Series Argument
    print("--- Step 3: The Taylor Series of the Fourier Transform ---")
    print("The Taylor series of F(xi) around xi=0 is given by the formula:")
    print("F(xi) = sum_{k=0 to inf} [F^(k)(0) / k!] * xi^k")
    print("Since every term F^(k)(0) is 0, the Taylor series is identically zero.")
    print("A key property of the Fourier transform of a Schwartz function is that it is an analytic function.")
    print("For an analytic function, if its Taylor series is zero, the function itself must be zero everywhere.")
    print("Therefore, F(xi) = 0 for all xi.\n")

    # Step 4: Conclusion from the Inverse Fourier Transform
    print("--- Step 4: Conclusion from the Inverse Fourier Transform ---")
    print("The function f(x) can be recovered from its Fourier transform F(xi) via the inverse transform:")
    print("f(x) = Integral(F(xi) * exp(2*pi*i*x*xi), (xi, -oo, oo))")
    print("Since we have concluded that F(xi) = 0, the integral above is zero.")
    print("This means f(x) = 0 for all x, which proves the statement.")

if __name__ == '__main__':
    solve_schwartz_moment_problem()