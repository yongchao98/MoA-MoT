import numpy as np
from scipy.integrate import quad

# This script illustrates the solution by demonstrating the contrapositive:
# If a Schwartz function f(x) is NOT the zero function, then at least one
# of its moments must be non-zero.

# We use the Gaussian function f(x) = exp(-x^2) as an example of a
# well-known, non-zero Schwartz function.
def f(x):
    """A non-zero Schwartz function (Gaussian)."""
    return np.exp(-x**2)

# The k-th moment is defined as M_k = integral from -inf to +inf of x^k * f(x) dx.
# We will now numerically compute the first few moments of our function.
print("Calculating moments for the non-zero Schwartz function f(x) = exp(-x^2):")
print("-" * 75)

# The exact value for the k=2n moment is Gamma(n + 1/2). For k=0, M_0 = sqrt(pi) ~= 1.772.
# We expect odd-numbered moments to be 0 and even-numbered moments to be non-zero.

for k in range(7):
    # This is the function we need to integrate to find the k-th moment.
    def integrand(x):
        return x**k * f(x)

    # Perform the numerical integration from -infinity to +infinity.
    # The quad function returns a tuple: (result, estimated_error).
    moment, error = quad(integrand, -np.inf, np.inf)

    # The problem asks to output each number in the final equation.
    # The equation is: Integral(x^k * f(x) dx) = moment.
    print(f"For k = {k}:")
    # Using an f-string to format the final equation with the calculated value.
    print(f"The equation is: integral from -inf to inf of (x^{k} * exp(-x^2)) dx = {moment:.8f}")

    if abs(moment) > 1e-9:
        print("Result: The moment is NON-ZERO.\n")
    else:
        # Treat values very close to zero as zero for clarity.
        print("Result: The moment is ZERO (within numerical precision).\n")

print("-" * 75)
print("As the calculations show, the non-zero function f(x) = exp(-x^2) has non-zero moments.")
print("This supports the conclusion that if all moments were zero, the function would have to be zero.")