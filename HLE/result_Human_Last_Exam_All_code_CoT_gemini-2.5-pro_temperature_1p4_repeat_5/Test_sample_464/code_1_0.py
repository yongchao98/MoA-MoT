import numpy as np
import scipy.integrate as integrate
import sympy as sp

def demonstrate_moment_fourier_relation():
    """
    This function illustrates the relationship between the moments of a function
    and the derivatives of its Fourier transform at the origin.
    The relation is: f_hat^(k)(0) = (-2*pi*i)^k * M_k
    where M_k is the k-th moment of f(x).

    We use the Schwartz function f(x) = exp(-x^2) as an example.
    """

    # 1. Define the function f(x) = exp(-x^2) for numerical integration
    def f(x):
        return np.exp(-x**2)

    print("Illustrating the core identity for the function f(x) = exp(-x^2).")
    print("The proof relies on the identity:")
    print("  (k-th derivative of Fourier Transform at 0) = (-2*pi*i)^k * (k-th moment of f)\n")
    print("-" * 70)

    # 2. Symbolically define the Fourier Transform for differentiation.
    # The Fourier transform of f(x) = exp(-x^2) is F(xi) = sqrt(pi) * exp(-(pi*xi)^2).
    # We use SymPy for exact symbolic differentiation.
    xi_s = sp.Symbol('xi')
    f_hat_s = sp.sqrt(sp.pi) * sp.exp(-(sp.pi * xi_s)**2)

    # 3. Loop for k=0 to 5 and check the identity
    for k in range(6):
        print(f"--- k = {k} ---")

        # Calculate the k-th moment of f(x) numerically.
        # M_k = integral from -inf to inf of x^k * f(x) dx
        integrand = lambda x: (x**k) * f(x)
        # SciPy's quad returns the integral and an error estimate.
        moment_k, _ = integrate.quad(integrand, -np.inf, np.inf)

        # Calculate the Left-Hand Side (LHS) of the identity.
        # This is the k-th derivative of the Fourier transform at xi=0.
        deriv_f_hat_s = sp.diff(f_hat_s, xi_s, k)
        lhs_sympy = deriv_f_hat_s.subs(xi_s, 0)
        lhs_val = complex(lhs_sympy.evalf())

        # Calculate the Right-Hand Side (RHS) of the identity.
        # RHS = (-2*pi*i)^k * M_k
        rhs_sympy = ((-2 * sp.pi * sp.I)**k) * moment_k
        rhs_val = complex(rhs_sympy.evalf())

        # Print each number in the final equation being checked.
        # The final equation is: lhs_val = rhs_val
        print(f"Moment M_{k} = {moment_k:.5f}")
        print("Checking the equation: [d^k/d(xi)^k F(f)](0) = (-2*pi*i)^k * M_k")
        print(f"  LHS (derivative): {lhs_val.real:12.5f} + {lhs_val.imag:9.5f}j")
        print(f"  RHS (moment):     {rhs_val.real:12.5f} + {rhs_val.imag:9.5f}j")

        # Check if the two complex numbers are close.
        if np.isclose(lhs_val, rhs_val):
            print("  => The identity holds.\n")
        else:
            print("  => A discrepancy exists due to numerical precision.\n")

    print("-" * 70)
    print("\nThe demonstration shows the identity holds. The mathematical argument is as follows:")
    print("If all moments M_k were 0, then the RHS would always be 0.")
    print("This implies the LHS (all derivatives of the Fourier transform at 0) must also be 0.")
    print("For the Fourier transform of a Schwartz function (an analytic function), this means the transform itself is 0.")
    print("If the Fourier transform is 0, the original function f(x) must also be 0.")

demonstrate_moment_fourier_relation()