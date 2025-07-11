import sympy as sp

def demonstrate_moment_problem_solution():
    """
    This script demonstrates the relationship between the moments of a function
    and the derivatives of its Fourier transform at the origin.
    """
    # Define symbolic variables
    x, xi = sp.symbols('x xi', real=True)
    f = sp.Function('f')(x)

    # Define the Fourier Transform using its integral definition.
    # The convention used is F(xi) = Integral(f(x)*exp(-2*pi*I*xi*x), (x, -oo, oo))
    fourier_transform_expr = sp.Integral(f * sp.exp(-2 * sp.I * sp.pi * xi * x), (x, -sp.oo, sp.oo))

    print("Let f(x) be a Schwartz function.")
    print("Let F(xi) be its Fourier transform, defined as:")
    print(f"F(xi) = {fourier_transform_expr}")
    print("\nLet's compute the derivatives of F(xi) with respect to xi, and evaluate them at xi = 0.")
    print("-" * 70)

    # Loop to compute derivatives for k = 0, 1, 2, 3
    for k in range(4):
        # Differentiate under the integral sign
        # sympy is smart enough to do this correctly
        derivative = sp.diff(fourier_transform_expr, xi, k)

        # Evaluate the derivative at xi = 0
        derivative_at_zero = derivative.subs(xi, 0)
        
        # Define the k-th moment
        moment_k = sp.Integral(x**k * f, (x, -sp.oo, sp.oo))
        
        print(f"For k = {k}:")
        
        # Show the k-th derivative of F(xi) evaluated at 0
        print(f"The {k}-th derivative of F(xi) at xi=0 is:")
        sp.pprint(derivative_at_zero)
        
        # Simplify the expression to show its relation to the k-th moment
        # Pull out the constant factor (-2*pi*I)**k from the integral
        constant_factor = (-2 * sp.I * sp.pi)**k
        simplified_expr = constant_factor * moment_k
        
        print("\nThis can be rewritten as:")
        equation = sp.Eq(sp.Symbol(f"F^({k})(0)"), simplified_expr)

        # sympy's pretty print can be a bit wide, let's format it ourselves.
        # equation for k=0: F^(0)(0) = Integral(f(x), (x, -oo, oo))
        # equation for k>0: F^(k)(0) = (-2*I*pi)**k * Integral(x**k*f(x), (x, -oo, oo))
        lhs = f"F^({k})(0)"
        rhs = f"{sp.pretty(constant_factor)} * {sp.pretty(moment_k)}"
        if k == 0:
            rhs = f"{sp.pretty(moment_k)}"

        print(f"{lhs} = {rhs}")
        print("\nIf the k-th moment Integral(x**k * f(x) dx) is 0, then F^({k})(0) = 0.")
        print("-" * 70)

    print("Since all moments are zero for k=0,1,2,..., all derivatives F^(k)(0) are zero.")
    print("Because F(xi) is analytic, its Taylor series around xi=0 is identically zero.")
    print("This implies F(xi) = 0 for all xi.")
    print("By the injectivity of the Fourier transform on Schwartz space, f(x) must also be 0.")

if __name__ == "__main__":
    demonstrate_moment_problem_solution()