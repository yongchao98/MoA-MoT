import sympy
from sympy import exp, pi, oo, I, sqrt

def demonstrate_moment_transform_identity():
    """
    Demonstrates the identity relating moments of a function to the
    derivatives of its Fourier transform at the origin.
    """
    # Set up symbolic variables for the function and its transform
    x, k_var = sympy.symbols('x k_var')

    # Define a non-zero Schwartz function, f(x) = exp(-x^2)
    f = exp(-x**2)
    print(f"Analyzing the function f(x) = {f}\n")

    # The Fourier Transform F(k) = integral(f*exp(-2*pi*i*x*k) dx)
    # For f(x)=exp(-x^2), its FT is sqrt(pi)*exp(-pi^2*k^2)
    f_hat = sympy.fourier_transform(f, x, k_var)
    print(f"Its Fourier Transform is F(k) = {f_hat}\n")

    print("Verifying the identity: F^(n)(0) = (-2*pi*i)^n * Moment_n(f)")
    print("-" * 70)

    # Check for the first few non-negative integers n
    for n in range(6):
        # 1. Calculate the n-th moment of f(x)
        # M_n = integral(x^n * f(x) dx) from -oo to oo
        moment_n = sympy.integrate(x**n * f, (x, -oo, oo))

        # 2. Calculate the Left-Hand Side (LHS) of the identity
        # The n-th derivative of the Fourier Transform F(k) w.r.t k, evaluated at k=0
        f_hat_deriv_n = sympy.diff(f_hat, k_var, n)
        lhs = f_hat_deriv_n.subs(k_var, 0)

        # 3. Calculate the Right-Hand Side (RHS) of the identity
        # (-2*pi*i)^n * M_n
        rhs = ((-2 * pi * I)**n * moment_n).simplify()

        # 4. Print the final equation with all values substituted
        print(f"For n = {n}:")
        print(f"Moment M_{n} = {moment_n}")
        print(f"FT derivative F^({n})(0) = {lhs}")
        # The user requested to print the final equation with numbers
        print(f"Resulting equation: {lhs} = {rhs}")
        print(f"Are sides equal? {sympy.simplify(lhs - rhs) == 0}\n")

    print("Conclusion:")
    print("The code demonstrates that if all moments M_n were 0 (as in the problem),")
    print("then all derivatives F^(n)(0) would also be 0. For an analytic function like")
    print("the Fourier transform, this implies F(k) is identically zero, which in turn")
    print("implies the original function f(x) is identically zero.")

if __name__ == '__main__':
    demonstrate_moment_transform_identity()
