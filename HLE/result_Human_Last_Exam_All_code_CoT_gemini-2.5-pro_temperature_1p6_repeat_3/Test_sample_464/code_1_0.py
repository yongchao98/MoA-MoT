import numpy as np
from scipy.integrate import quad
from scipy.misc import derivative

def run_demonstration():
    """
    Demonstrates the relationship between moments of a function and the
    derivatives of its Fourier transform.
    """
    print("This script demonstrates the proof that if all moments of a Schwartz function are zero,")
    print("the function itself must be zero.\n")
    print("The key equation is: D_k = (-2 * pi * i)^k * M_k")
    print("Where M_k is the k-th moment of f(x) and D_k is the k-th derivative of its Fourier transform at 0.\n")
    print("If all moments M_k are 0, then all derivatives D_k must be 0.")
    print("This implies the Fourier transform is the zero function, and therefore f(x) is the zero function.\n")

    # We use a non-zero Schwartz function, f(x) = exp(-x^2), to show the relationship holds.
    # Its moments are NOT all zero.

    # The Schwartz function f(x) = exp(-x^2)
    def f(x):
        return np.exp(-x**2)

    # Analytically known Fourier transform of f(x) is F(xi) = sqrt(pi) * exp(-(pi*xi)^2)
    def f_hat(xi):
        return np.sqrt(np.pi) * np.exp(-(np.pi * xi)**2)

    print("Demonstration for the non-zero function f(x) = exp(-x^2):")
    print("-" * 75)
    print("{:<3} | {:<20} | {:<25} | {:<25}".format("k", "Moment M_k", "FT Deriv D_k at 0", "Check: (-2*pi*i)^k*M_k"))
    print("-" * 75)

    for k in range(5):
        # Calculate the k-th moment: M_k = integral(x^k * f(x) dx)
        moment_integrand = lambda x: (x**k) * f(x)
        # The integral of an odd function over a symmetric interval is 0.
        # quad might have small numerical errors.
        if k % 2 != 0:
            moment_val = 0.0
        else:
            moment_val, _ = quad(moment_integrand, -np.inf, np.inf)

        # Calculate the k-th derivative of the Fourier transform at xi=0
        # For an even/odd function, the derivatives of its real-valued FT will be zero or non-zero.
        if k % 2 != 0:
            ft_deriv_val = 0.0
        else:
            # We must use an odd order for the central difference formula of derivative, and order > n.
            order = k + 1 + (k % 2)
            ft_deriv_val = derivative(f_hat, 0.0, dx=1e-5, n=k, order=order)

        # Check the relationship: D_k = (-2*pi*i)^k * M_k
        check_val = ((-2j * np.pi)**k) * moment_val
        
        # We check the real part, as the FT is real-valued, so its derivatives are real.
        print("{:<3} | {:<20.6f} | {:<25.6f} | {:<25.6f}".format(k, moment_val, ft_deriv_val, check_val.real))
    
    print("-" * 75)

if __name__ == '__main__':
    run_demonstration()
