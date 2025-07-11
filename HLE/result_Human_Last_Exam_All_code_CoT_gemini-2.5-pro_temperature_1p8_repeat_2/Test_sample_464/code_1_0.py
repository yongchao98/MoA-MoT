import numpy as np
from scipy.integrate import quad
from scipy.misc import derivative

def demonstrate_moment_transform_identity():
    """
    This function demonstrates the mathematical identity linking a function's
    moments to the derivatives of its Fourier transform, which is central
    to proving that a Schwartz function with all-zero moments must be zero.

    We use the Schwartz function f(x) = exp(-pi*x^2). Its Fourier transform is
    f_hat(xi) = exp(-pi*xi^2).
    The identity is: f_hat^(k)(0) = (-2*pi*i)^k * Integral(x^k * f(x) dx)
    """

    # Define the Schwartz function f(x) = exp(-pi*x^2)
    f = lambda x: np.exp(-np.pi * x**2)

    # The analytical Fourier Transform of f(x) is f_hat(xi) = exp(-pi*xi^2)
    f_hat = lambda xi: np.exp(-np.pi * xi**2)

    print("Demonstrating the identity: (d^k/dxi^k)f_hat(0) = (-2*pi*i)^k * moment_k\n")
    print("-" * 70)

    # We will check the identity for k = 0, 1, 2, 3
    for k in range(4):
        print(f"Case k = {k}:")

        # 1. Calculate the k-th moment: M_k = Integral(x^k * f(x) dx)
        integrand = lambda x: (x**k) * f(x)
        # quad returns a tuple (result, error_bound)
        moment_k, _ = quad(integrand, -np.inf, np.inf)

        # 2. Calculate the k-th derivative of the Fourier transform at xi=0
        # For complex functions, derivative works by looking at the real part.
        # Since f_hat is real, this is fine.
        deriv_f_hat_at_0 = derivative(f_hat, 0.0, dx=1e-6, n=k, order=k+1+(k%2))
        
        # 3. Calculate the right-hand side (RHS) of the identity
        rhs = ((-2j * np.pi)**k) * moment_k

        # 4. Print the numbers for comparison
        print(f"  k-th moment (integral): M_{k} = {moment_k:.8f}")
        print(f"  LHS: k-th derivative of f_hat at 0 = {deriv_f_hat_at_0:.8f}")
        # Note: We take the real part of RHS since LHS is real for our choice of f_hat.
        # In the general case, both can be complex.
        print(f"  RHS: (-2*pi*i)^{k} * M_{k} = {rhs.real:.8f} + {rhs.imag:.8f}i")

        # Conclusion for this k
        if np.isclose(deriv_f_hat_at_0, rhs.real) and np.isclose(0, rhs.imag):
            print("  LHS and RHS match.")
            # For the hypothetical case where f has all zero moments
            if k == 0:
                print("  If all moments were 0, M_k would be 0, forcing the derivative to be 0.")
        else:
            print("  LHS and RHS DO NOT match.")
        
        print("-" * 70)
        
demonstrate_moment_transform_identity()