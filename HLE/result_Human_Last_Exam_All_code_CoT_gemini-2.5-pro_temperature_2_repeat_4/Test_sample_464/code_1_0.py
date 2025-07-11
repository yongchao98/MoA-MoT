import numpy as np
from scipy.integrate import quad
from scipy.special import gamma

def schwartz_function_illustration():
    """
    Calculates the moments of a non-zero Schwartz function, f(x) = exp(-x^2),
    to demonstrate that its moments are not all zero.
    """

    # Define a non-zero Schwartz function, the Gaussian function.
    f = lambda x: np.exp(-x**2)

    print("Illustrating the contrapositive: a non-zero Schwartz function f(x) = exp(-x^2) has non-zero moments.")
    print("A moment is defined by the equation: M_k = integral from -inf to inf of x^k * f(x) dx.")
    print("-" * 70)

    # Calculate moments for k from 0 to 10
    for k in range(11):
        # Define the function to integrate for the k-th moment
        integrand = lambda x: (x**k) * f(x)

        # Perform numerical integration from -infinity to +infinity
        moment_val, error = quad(integrand, -np.inf, np.inf)

        # For f(x) = exp(-x^2), moments for odd k are 0, and for even k=2n are Gamma(n + 1/2).
        if k % 2 == 1:
            analytical_val = 0.0
        else:
            n = k / 2.0
            analytical_val = gamma(n + 0.5)

        # Print the final equation with the calculated numbers
        # This shows the equation form: Integral = result
        print(f"For k={k}:")
        print(f"  The equation is: integral(x^{k} * exp(-x^2)) dx = {moment_val:.8f}")
        print(f"  The known analytical value is: {analytical_val:.8f}")

    print("-" * 70)
    print("As the output shows, the moments for even values of k are non-zero.")
    print("This is consistent with the theorem: since the function is not identically zero, its moments are not all zero.")

schwartz_function_illustration()